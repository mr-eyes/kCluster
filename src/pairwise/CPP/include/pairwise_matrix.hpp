#include <vector>
#include <stdint.h>

class Matrix{
private:
    uint32_t *** matrix = NULL;
    int SIZE_XY;
    int SIZE_Z;

public:
    Matrix(int tr_no, int Qs) {
        SIZE_XY = tr_no;
        SIZE_Z = Qs;
        matrix = new uint32_t **[tr_no];
        for (int i = 0; i < tr_no; i++) {
            matrix[i] = new uint32_t *[tr_no];
            for (int j = 0; j < tr_no; j++) {
                matrix[i][j] = new uint32_t[Qs];
                for (int k = 0; k < Qs; k++) {
                    matrix[i][j][k] = 0;
                }
            }
        }
    }

    void set(uint32_t seq1, uint32_t seq2, uint8_t Q, uint32_t val){
        matrix[seq1-1][seq2-1][Q-1] += val;
    }

    uint32_t get(uint32_t seq1, uint32_t seq2, uint8_t Q){
        return matrix[seq1-1][seq2-1][Q-1];
    }

    void dump(){
        for(int Seq1 = 0; Seq1 < SIZE_XY; ++Seq1)
        {
            for(int Seq2 =0; Seq2 < SIZE_XY; ++Seq2)
            {
                for(int Q = 0; Q < SIZE_Z; ++Q)
                {
                    if (matrix[Seq1][Seq2][Q])
                    cout << "pair["<< Seq1 + 1 << "][" << Seq2 + 1 << "][" << Q + 1 << "] = " << matrix[Seq1][Seq2][Q] << endl;
                }
            }
        }
    }

};