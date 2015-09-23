#include <iostream>
#include <cstdlib>
#include <string>

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

int * classnumbers;
long classnumber_tablesize;

int init_classnumbers() {
    std::string home = std::getenv("HOME");
    std::string filelocation = home + "/include/classnumbers";
    int fd = open(filelocation.c_str(), O_RDONLY);
    if(fd == -1) {
        std::cerr << "error opening mmapped file" << std::endl;
        exit(-1);
    }
    struct stat fs;
    if(fstat(fd, &fs) == -1) {
        std::cerr << "error in fstat" << std::endl;
        exit(-1);
    }

    size_t filesize = fs.st_size;
    classnumber_tablesize = filesize/sizeof(int);

    classnumbers = (int *)mmap(NULL, filesize, PROT_READ, MAP_SHARED, fd, 0);
    if(classnumbers == MAP_FAILED) {
        std::cerr << "error opening mmapped file" << std::endl;
        exit(-1);
    }
    return 1;
}
