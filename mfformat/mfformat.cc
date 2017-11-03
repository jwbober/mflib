#include <string>

#include <cstdio>
#include <cstdint>
#include <cstring>

#include "mfformat.h"
#include "flint/fmpz_poly.h"
#include "slint.h"

#include "sqlite3.h"

extern "C" {

using namespace std;

int write_mfheader(FILE * outfile, struct mfheader * header) {
    if(!fwrite((char*)&header->version, sizeof(header->version), 1, outfile)) return 0;
    if(!fwrite((char*)&header->level, sizeof(header->level), 1, outfile)) return 0;
    if(!fwrite((char*)&header->weight, sizeof(header->weight), 1, outfile)) return 0;
    if(!fwrite((char*)&header->chi, sizeof(header->chi), 1, outfile)) return 0;
    if(!fwrite((char*)&header->orbit, sizeof(header->orbit), 1, outfile)) return 0;
    if(!fwrite((char*)&header->j, sizeof(header->j), 1, outfile)) return 0;
    if(!fwrite((char*)&header->prec, sizeof(header->prec), 1, outfile)) return 0;
    if(!fwrite((char*)&header->exponent, sizeof(header->exponent), 1, outfile)) return 0;
    if(!fwrite((char*)&header->ncoeffs, sizeof(header->ncoeffs), 1, outfile)) return 0;
    if(!fwrite((char*)&header->reserved, sizeof(header->reserved), 1, outfile)) return 0;
    return 1;
}

int read_mfheader(FILE * outfile, struct mfheader * header) {
    if(!fread((char*)&header->version, sizeof(header->version), 1, outfile)) return 0;
    if((header->version != MFV1) && (header->version != MFV2)) return 0;
    if(!fread((char*)&header->level, sizeof(header->level), 1, outfile)) return 0;
    if(!fread((char*)&header->weight, sizeof(header->weight), 1, outfile)) return 0;
    if(!fread((char*)&header->chi, sizeof(header->chi), 1, outfile)) return 0;
    if(!fread((char*)&header->orbit, sizeof(header->orbit), 1, outfile)) return 0;
    if(!fread((char*)&header->j, sizeof(header->j), 1, outfile)) return 0;
    if(!fread((char*)&header->prec, sizeof(header->prec), 1, outfile)) return 0;
    if(!fread((char*)&header->exponent, sizeof(header->exponent), 1, outfile)) return 0;
    if(!fread((char*)&header->ncoeffs, sizeof(header->ncoeffs), 1, outfile)) return 0;
    if(!fread((char*)&header->reserved, sizeof(header->reserved), 1, outfile)) return 0;
    return 1;
}

void acb_set_mfcoeff(acb_t out, fmpz_t x, fmpz_t y, struct mfheader * header) {
    acb_set_fmpz_fmpz(out, x, y);
    acb_mul_2exp_si(out, out, header->exponent);
    if(header->prec != MF_PREC_EXACT) {
        mag_set_ui(arb_radref(acb_realref(out)), 1);
        mag_mul_2exp_si(arb_radref(acb_realref(out)), arb_radref(acb_realref(out)), header->prec);
        if(header->chi != 1) {
            mag_set_ui(arb_radref(acb_imagref(out)), 1);
            mag_mul_2exp_si(arb_radref(acb_imagref(out)), arb_radref(acb_imagref(out)), header->prec);
        }
    }
}

size_t acb_write_mfcoeff(FILE * outfile, struct mfheader * header, acb_t coeff) {
    fmpz_t x;
    fmpz_t y;
    acb_t z;

    fmpz_init(x);
    fmpz_init(y);
    acb_init(z);

    size_t retval = 0;

    acb_mul_2exp_si(z, coeff, -header->exponent);
    if( !arb_get_unique_fmpz(x, acb_realref(z)) ) {
        arb_floor(acb_realref(z), acb_realref(z), -header->prec + 500);
        if( !arb_get_unique_fmpz(x, acb_realref(z)) ) {
            goto cleanup;
        }
    }

    if( !arb_get_unique_fmpz(y, acb_imagref(z)) ) {
        arb_floor(acb_imagref(z), acb_imagref(z), -header->prec + 500);
        if( !arb_get_unique_fmpz(y, acb_imagref(z)) ) {
            goto cleanup;
        }
    }

    retval += fmpz_out_raw(outfile, x);
    if(header->chi != 1)
        retval += fmpz_out_raw(outfile, y);

cleanup:
    fmpz_clear(x);
    fmpz_clear(y);
    acb_clear(z);

    return retval;
}

int acb_attempt_rounding_mfcoeff(struct mfheader * header, acb_t coeff) {
    fmpz_t x;
    fmpz_t y;
    acb_t z;

    fmpz_init(x);
    fmpz_init(y);
    acb_init(z);

    int retval = 1;

    acb_mul_2exp_si(z, coeff, -header->exponent);
    if( !arb_get_unique_fmpz(x, acb_realref(z)) ) {
        arb_floor(acb_realref(z), acb_realref(z), -header->prec + 500);
        if( !arb_get_unique_fmpz(x, acb_realref(z)) ) {
            retval = 0;
            goto cleanup;
        }
    }

    if( !arb_get_unique_fmpz(y, acb_imagref(z)) ) {
        arb_floor(acb_imagref(z), acb_imagref(z), -header->prec + 500);
        if( !arb_get_unique_fmpz(y, acb_imagref(z)) ) {
            retval = 0;
            goto cleanup;
        }
    }

cleanup:
    fmpz_clear(x);
    fmpz_clear(y);
    acb_clear(z);

    return retval;
}




int read_mffile(FILE * infile, struct mfheader * header, acb_ptr * coeffs) {
    //
    if(read_mfheader(infile, header) == 0) return 0;
    return read_mfdata(infile, header, coeffs);
}

int read_mfdata(FILE * infile, struct mfheader * header, acb_ptr * coeffs) {
    if(header->version == MFV2) {
        *coeffs = _acb_vec_init(header->ncoeffs);
        if(!*coeffs) return 0;
        for(int k = 0; k < header->ncoeffs; k++) {
            arb_set_str(acb_realref(*coeffs + k), "nan", 10);
        }

        fmpz_t x;
        fmpz_t y;
        fmpz_init(x);
        fmpz_init(y);

        int k = 0;
        int p = prime_powers_table[k];
        while(p <= header->ncoeffs) {
            fmpz_inp_raw(x, infile);
            if(header->chi != 1) {
                fmpz_inp_raw(y, infile);
            }
            acb_set_mfcoeff((*coeffs) + p - 1, x, y, header);
            k++;
            p = prime_powers_table[k];
        }
        acb_set_ui(*coeffs, 1);
        k = 0;
        int p1 = prime_powers_table[k];
        k = 0;
        while(p1 <= header->ncoeffs) {
            for(int j = 2; p1*j <= header->ncoeffs; j++) {
                if(GCD(p1, j) != 1)
                    continue;
                acb_mul(*coeffs + p1*j - 1, *coeffs + p1 - 1, *coeffs + j - 1, 1000);
            }
            k++;
            p1 = prime_powers_table[k];
        }
        return 1;
    }

    *coeffs = _acb_vec_init(header->ncoeffs);
    if(!*coeffs) return 0;

    fmpz_t x;
    fmpz_t y;
    fmpz_init(x);
    fmpz_init(y);

    for(unsigned int k = 0; k < header->ncoeffs; k++) {
        fmpz_inp_raw(x, infile);
        if(header->chi != 1) {
            fmpz_inp_raw(y, infile);
        }
        acb_set_mfcoeff((*coeffs) + k, x, y, header);
    }

    return 1;
}

int read_mfdatablob(const void * data, int datasize, struct mfheader * header, acb_ptr * coeffs) {
    FILE * f = fmemopen((void *)data, datasize, "r");
    read_mfdata(f, header, coeffs);
    fclose(f);
    return 1;
}

int write_mffile(FILE * outfile, struct mfheader * header, acb_ptr coeffs) {
    if(header->version != MFV2) return 0;
    if(!write_mfheader(outfile, header)) return 0;
    return write_mfdata(outfile, header, coeffs);
}

int write_mfdata(FILE * outfile, struct mfheader * header, acb_ptr coeffs) {
    if(header->version != MFV2) return 0;
    int k = 0;
    int p = prime_powers_table[k];
    while(p <= header->ncoeffs) {
        if(acb_write_mfcoeff(outfile, header, coeffs + p - 1) == 0) {
            return 0;
        }
        k++;
        p = prime_powers_table[k];
        if(p == 0) return 0;
    }
    return 1;
}

int insert_into_sqlite(sqlite3 * db, struct mfheader * header, acb_ptr coeffs) {
    // expect a database with a table modforms
    // CREATE TABLE modforms (level INTEGER, weight INTEGER, chi INTEGER, orbit INTEGER, j INTEGER,
    //      prec INTEGER, exponent INTEGER, ncoeffs INTEGER, coefficients BLOB)

    char * filedata;
    size_t filesize;
    FILE * coeff_file = open_memstream(&filedata, &filesize);
    fclose(coeff_file);

    sqlite3_stmt * stmt;
    char sql[] = "INSERT INTO modforms (level, weight, chi, orbit, j, prec, exponent, ncoeffs, coefficients) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);";

    int result = sqlite3_prepare_v2(db, sql, sizeof(sql), &stmt, NULL);
    if(result != SQLITE_OK) return result;
    sqlite3_bind_int(stmt, 1, header->level);
    sqlite3_bind_int(stmt, 2, header->weight);
    sqlite3_bind_int(stmt, 3, header->chi);
    sqlite3_bind_int(stmt, 4, header->orbit);
    sqlite3_bind_int(stmt, 5, header->j);
    sqlite3_bind_int(stmt, 6, header->prec);
    sqlite3_bind_int(stmt, 7, header->exponent);
    sqlite3_bind_int(stmt, 8, header->ncoeffs);
    sqlite3_bind_blob(stmt, 9, (void *)filedata, filesize, free);

    result = sqlite3_step(stmt);
    sqlite3_finalize(stmt);
    return result;
}

void iterate_through_sqlitefile(
        const char * filename,
        int (*callback)(struct mfheader * header, int coeff_datasize, const void * coeff_data, acb_ptr coeffs),
        int populate_coefficients) {
    sqlite3 * db;
    sqlite3_open(filename, &db);
    char sql[] = "SELECT level, weight, chi, orbit, j, prec, exponent, ncoeffs, coefficients FROM modforms";
    sqlite3_stmt * stmt;
    int result = sqlite3_prepare_v2(db, sql, sizeof(sql), &stmt, NULL);
    struct mfheader header;
    const void * coeff_data;
    int coeff_datasize;
    acb_ptr coeffs;
    int callback_result = 0;

    result = sqlite3_step(stmt);
    while(result == SQLITE_ROW && callback_result == 0) {
        header.level = sqlite3_column_int(stmt, 0);
        header.weight = sqlite3_column_int(stmt, 1);
        header.chi = sqlite3_column_int(stmt, 2);
        header.orbit = sqlite3_column_int(stmt, 3);
        header.j = sqlite3_column_int(stmt, 4);
        header.prec = sqlite3_column_int(stmt, 5);
        header.exponent = sqlite3_column_int(stmt, 6);
        header.ncoeffs = sqlite3_column_int(stmt, 7);
        coeff_data = sqlite3_column_blob(stmt, 8);
        coeff_datasize = sqlite3_column_bytes(stmt, 8);

        header.version = MFV2;

        if(populate_coefficients) {
            read_mfdatablob(coeff_data, coeff_datasize, &header, &coeffs);
        }

        callback_result = callback(&header, coeff_datasize, coeff_data, coeffs);

        if(populate_coefficients) {
            _acb_vec_clear(coeffs, header.ncoeffs);
        }
        if(callback_result == 0) result = sqlite3_step(stmt);
    }

    sqlite3_finalize(stmt);
    sqlite3_close(db);
}

void iterate_through_sqlitefile_with_filter(
        const char * filename,
        int (*callback)(struct mfheader * header, int coeff_datasize, const void * coeff_data, acb_ptr coeffs),
        int populate_coefficients,
        int level,
        int weight,
        int chi) {
    string where_clause;
    if(level == 0 && weight == 0 && chi == 0) {
        where_clause = "";
    }
    else {
        where_clause = " WHERE ";
    }
    bool insert_and = false;
    if(level != 0) {
        where_clause += " level=" + to_string(level) + " ";
        insert_and = true;
    }
    if(weight != 0) {
        if(insert_and) where_clause += " AND ";
        where_clause += " weight=" + to_string(weight) + " ";
        insert_and = true;
    }
    if(chi != 0) {
        if(insert_and) where_clause += " AND ";
        where_clause += " chi=" + to_string(chi) + " ";
    }

    sqlite3 * db;
    sqlite3_open(filename, &db);
    string sql = "SELECT level, weight, chi, orbit, j, prec, exponent, ncoeffs, coefficients FROM modforms " + where_clause;
    sqlite3_stmt * stmt;
    int result = sqlite3_prepare_v2(db, sql.c_str(), sql.size(), &stmt, NULL);
    struct mfheader header;
    const void * coeff_data;
    int coeff_datasize;
    acb_ptr coeffs;
    int callback_result = 0;

    result = sqlite3_step(stmt);
    while(result == SQLITE_ROW && callback_result == 0) {
        header.level = sqlite3_column_int(stmt, 0);
        header.weight = sqlite3_column_int(stmt, 1);
        header.chi = sqlite3_column_int(stmt, 2);
        header.orbit = sqlite3_column_int(stmt, 3);
        header.j = sqlite3_column_int(stmt, 4);
        header.prec = sqlite3_column_int(stmt, 5);
        header.exponent = sqlite3_column_int(stmt, 6);
        header.ncoeffs = sqlite3_column_int(stmt, 7);
        coeff_data = sqlite3_column_blob(stmt, 8);
        coeff_datasize = sqlite3_column_bytes(stmt, 8);

        header.version = MFV2;

        if(populate_coefficients) {
            read_mfdatablob(coeff_data, coeff_datasize, &header, &coeffs);
        }

        callback_result = callback(&header, coeff_datasize, coeff_data, coeffs);

        if(populate_coefficients) {
            _acb_vec_clear(coeffs, header.ncoeffs);
        }
        if(callback_result == 0) result = sqlite3_step(stmt);
    }

    sqlite3_finalize(stmt);
    sqlite3_close(db);
}




int mfdb_get_entry(sqlite3 * db, struct mfheader * header, acb_ptr * coeffs, int level, int weight, int chi, int j) {
    char sql[] = "SELECT level, weight, chi, orbit, j, prec, exponent, ncoeffs, coefficients FROM modforms "
                    "WHERE level=? and weight=? and chi=? and j=? LIMIT 1";
    sqlite3_stmt * stmt;
    int result = sqlite3_prepare_v2(db, sql, sizeof(sql), &stmt, NULL);
    sqlite3_bind_int(stmt, 1, level);
    sqlite3_bind_int(stmt, 2, weight);
    sqlite3_bind_int(stmt, 3, chi);
    sqlite3_bind_int(stmt, 4, j);

    result = sqlite3_step(stmt);
    if(result != SQLITE_ROW) {
        return 0;
    }
    header->level = sqlite3_column_int(stmt, 0);
    header->weight = sqlite3_column_int(stmt, 1);
    header->chi = sqlite3_column_int(stmt, 2);
    header->orbit = sqlite3_column_int(stmt, 3);
    header->j = sqlite3_column_int(stmt, 4);
    header->prec = sqlite3_column_int(stmt, 5);
    header->exponent = sqlite3_column_int(stmt, 6);
    header->ncoeffs = sqlite3_column_int(stmt, 7);

    const void * coeff_data = sqlite3_column_blob(stmt, 8);
    size_t coeff_datasize = sqlite3_column_bytes(stmt, 8);

    header->version = MFV2;
    read_mfdatablob(coeff_data, coeff_datasize, header, coeffs);

    sqlite3_finalize(stmt);
    return 1;
}

} // extern "C"
