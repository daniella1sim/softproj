#ifndef SYMNMF_H_
#define SYMNMF_H_

struct cord
{
    double value;
    struct cord *next;
};

struct vector
{
    struct vector *next;
    struct cord *cords;
};

typedef struct
{
    double **data;
    int rows;
    int cols;
} Matrix;


Matrix *initializeMatrix(int n, int m);
Matrix *transpose(Matrix * matrix);
double MatrixDistance(Matrix *matrixA, Matrix *matrixB);
Matrix *matrixMultiply(Matrix *matrix1, Matrix *matrix2);

struct cord* createNewCord();
struct vector* createNewVector();
struct cord* addNewCord(struct cord *currCord);
struct vector* addVector(struct vector *currVec);
void freeVector(struct vector *v);
int countVectors(struct vector *headVec);

void freeMatrix(Matrix *matrix);
void printMatrix(Matrix *matrix);

Matrix *similarityMatrix(struct vector *points, int numOfPoints);
Matrix *diagonalDegreeMatrix(Matrix *similarityMatrix);
Matrix *normalizedSimilarityMatrix(Matrix *similarityMat, Matrix *diagonalDegreeMat);

#endif
