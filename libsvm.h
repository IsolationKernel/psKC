#ifndef LIBSVM_H
#define LIBSVM_H

class csr_c;
class file_name_c;

class libsvm_c {
    public:
        libsvm_c(file_name_c const & file_name);

        void load(csr_c & data);

    private:
        file_name_c const & file_name;
};

#endif // LIBSVM_H
