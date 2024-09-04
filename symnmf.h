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

struct cord* createNewCord();
struct vector* createNewVector();
struct cord* addNewCord(struct cord *currCord);
struct vector* addVector(struct vector *currVec);
void freeVector(struct vector *v);
int countVectors(struct vector *headVec);

void freeMatrix(double **matrix, int n);
void printMatrix(double** matrix, int n);

double** similarityMatrix(struct vector *points, int numOfPoints);
double** diagonalDegreeMatrix(double **similarityMatrix, int n);
double** normalizedSimilarityMatrix(double **similarityMat, double **diagonalDegreeMat, int n);

#endif