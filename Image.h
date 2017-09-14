#ifndef PARALLEL_STRATEGIES_TESTING_IMAGE_H
#define PARALLEL_STRATEGIES_TESTING_IMAGE_H
#include <fstream> // for file I/O

using namespace std;
typedef unsigned char unchar; // Easier to understand & code.

class MImage {

public:
    MImage(const char* fileName);
    ~MImage(); //Deconstructor
    void write(const char* fileName);
    unchar** getData() const;
private:
    ifstream* pInFile;
    ofstream* pOutFile;
    unchar imageHeaderData[1078];
    unchar** imageData;
    unchar m_smoothFilter[3][3];
    unchar**  filteredData;
};

#endif //PARALLEL_STRATEGIES_TESTING_IMAGE_H
