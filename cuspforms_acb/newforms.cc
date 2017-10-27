#include <iostream>
#include <vector>
#include <ctime>

#include "cuspforms_acb.h"
#include "arb-extras.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::clock;

void cuspforms_acb::newforms(acb_mat_t out, int ncoeffs) {
    int dim = new_dimension();
    if(verbose) {
        cout << clock()/CLOCKS_PER_SEC << ": computing newform basis for space of dimension " << dim << endl;
    }
    if(dim == 0) {acb_mat_init(out, 0, 0);}
    acb_t z1;
    acb_init(z1);
    vector<int> basis_data = newspace_basis_data();
    if(verbose > 1) {
        cout << clock()/CLOCKS_PER_SEC << ": basis n for newspace:";
        for(int n : basis_data) cout << " " << n;
        cout << endl;
    }
    int coefficients_needed_for_full_rank = basis_data[dim - 1];
    ncoeffs = std::max(ncoeffs, 2 * coefficients_needed_for_full_rank + 5);
    if(verbose) cout << clock()/CLOCKS_PER_SEC << ": computing traces" << endl;
    compute_traces(ncoeffs * basis_data[dim-1]);
    if(verbose) cout << clock()/CLOCKS_PER_SEC << ": finished initial computation of traces" << endl;

    int sqrt_lim = std::max(1000, basis_data[dim - 1] + 1); // XXX. This is a potential problem.
                                                            // There is no reason that p in the
                                                            // !found_unique_eigenvalues loop could
                                                            // not get much bigger. (But N will have
                                                            // go get really large for that to happen,
                                                            // probably.
    acb_ptr invsqrts = _acb_vec_init(sqrt_lim);
    acb_set_ui(&invsqrts[1], 1u);
    for(int p = 2; p < sqrt_lim; p = next_prime(p)) {
        acb_reasonable_sqrt(&invsqrts[p], chi_values[p % level], prec);
        acb_conj(&invsqrts[p], &invsqrts[p]);
        for(int k = p; k < sqrt_lim; k += p) {
            acb_mul(&invsqrts[k], &invsqrts[k/p], &invsqrts[p], prec);
        }
    }


    arb_mat_t basis; // Will actually be a subset of the fourier coefficients
                     // of a basis. This will be a full rank symmetric matrix
                     // with entry (i,j) equal to T_{n_i, n_j}, where
                     // n_i is the ith entry in basis_data

    arb_mat_init(basis, dim, dim);
    for(int k = 0; k < dim; k++) {
        for(int j = k; j < dim; j++) {
            int n = basis_data[k];
            int m = basis_data[j];
            trace_TnTm(z1, n, m);
            acb_mul(z1, z1, &invsqrts[n], prec);
            acb_mul(z1, z1, &invsqrts[m], prec);
            arb_set(arb_mat_entry(basis, j, k), acb_realref(z1));
            arb_set(arb_mat_entry(basis, k, j), acb_realref(z1));
        }
    }


    if(verbose) cout << clock()/CLOCKS_PER_SEC << ": computed basis matrix" << endl;

    arb_mat_t Tp_basis; // The basis after being acted on by Tp (or a sum of Tp)
    arb_mat_init(Tp_basis, dim, dim);
    int p = 2;
    bool found_unique_eigenvalues = false;

    acb_mat_t basis_transformation;
    arb_mat_t eigenvectors;
    arb_mat_t eigenvalues;

    acb_mat_init(basis_transformation, dim, dim);
    arb_mat_init(eigenvalues, dim, 1);
    arb_mat_init(eigenvectors, dim, dim);

    //rmatrix_t rr_Tp_basis(dim, dim);
    bool first = true;
    while(!found_unique_eigenvalues) {
        while(GCD(p, level) > 1) {p++;}
        if(verbose) cout << clock()/CLOCKS_PER_SEC << ": computing hecke action for p = " << p << endl;
        if(verbose) cout << clock()/CLOCKS_PER_SEC << ": (maybe) computing more traces." << endl;
        compute_traces((p+1) * (basis_data[dim-1] + 1)*(basis_data[dim-1] + 1) + 5);

        if(verbose) cout << clock()/CLOCKS_PER_SEC << ": preparing rotated basis matrix." <<  endl;

        int multiplier = 1;
        if(!first) {
            multiplier = ((double)rand()/(double)RAND_MAX) * 5 + 1;
        }
        first = false;

        for(int k = 0; k < dim; k++) {
            for(int j = k; j < dim; j++) {
                int n = basis_data[k];
                int m = basis_data[j];
                trace_TpTnTm(z1, p, n, m);
                acb_mul(z1, z1, &invsqrts[n], prec);
                acb_mul(z1, z1, &invsqrts[m], prec);
                acb_mul(z1, z1, &invsqrts[p], prec);
                //arb_add(arb_mat_entry(Tp_basis, j, k), arb_mat_entry(Tp_basis, j, k), acb_realref(z1), prec);
                arb_addmul_si(arb_mat_entry(Tp_basis, j, k), acb_realref(z1), multiplier, prec);
                if(j != k) {
                    //arb_add(arb_mat_entry(Tp_basis, k, j), arb_mat_entry(Tp_basis, k, j), acb_realref(z1), prec);
                    arb_addmul_si(arb_mat_entry(Tp_basis, k, j), acb_realref(z1), multiplier, prec);
                }
            }
        }

        if(verbose > 1)
            arb_mat_printd(Tp_basis, 10);
        if(verbose > 1)
            cout << endl;

        if(verbose) cout << clock()/CLOCKS_PER_SEC << ": solving eigenvalue problem." <<  endl;
        int result = arb_mat_generalized_eigenproblem_symmetric_positive_definite(eigenvalues, eigenvectors, Tp_basis, basis, prec);
        if(verbose) cout << clock()/CLOCKS_PER_SEC << ": finished eigenvalue problem." <<  endl;

        if(verbose > 1) {
            cout << clock()/CLOCKS_PER_SEC << ": eigenvalues:" << endl;
            arb_mat_print_sage_float(eigenvalues);
            cout << endl;
            arb_mat_printd(eigenvalues, prec/3);
            cout << endl;
            cout << clock()/CLOCKS_PER_SEC << ": eigenvectors:" << endl;
            arb_mat_print_sage_float(eigenvectors);
            cout << endl;
            arb_mat_printd(eigenvectors, 10);
            cout << endl;
        }

        if(result == 0) {
            for(int j = 0; j < dim; j++) {
                for(int k = 0; k < dim; k++) {
                    acb_mul_arb(    acb_mat_entry(basis_transformation, j, k),
                                    &invsqrts[basis_data[j]],
                                    arb_mat_entry(eigenvectors, j, k),
                                    prec);
                }
            }
            found_unique_eigenvalues = true;
        }
        else {
            p++;
            if(verbose) {
                cout << "using more Hecke operators to get unique eigenvalues. next p = " << p << endl;
            }
        }

    }

    if(verbose) cout << clock()/CLOCKS_PER_SEC << ": finished computation of traces." << endl;

    newspace_basis(out, 2);
    acb_mat_mul(out, out, basis_transformation, prec);
    for(int j = 0; j < dim; j++) {
        for(int k = 0; k < dim; k++) {
            acb_div(acb_mat_entry(basis_transformation, k,j), acb_mat_entry(basis_transformation,k,j), acb_mat_entry(out,1, j), prec);
        }
    }
    acb_mat_clear(out);
    newspace_basis(out, ncoeffs);
    if(verbose) {
        cout << clock()/CLOCKS_PER_SEC << ": newforms_acb::newforms(): final matrix multiplication." << endl;
    }
    acb_mat_mul(out, out, basis_transformation, prec);
    if(verbose) {
        cout << clock()/CLOCKS_PER_SEC << ": newforms_acb::newforms(): finished final matrix multiplication." << endl;
    }
}
