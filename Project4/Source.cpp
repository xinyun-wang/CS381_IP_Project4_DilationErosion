#include <fstream>
#include <iostream>
using namespace std;

class morphology {
public:
    int numImgR;
    int numImgC;
    int imgMin;
    int imgMax;
    int numStructR;
    int numStructC;
    int structMin;
    int structMax;
    int rowOri;
    int colOri;
    int rowFrameSize;
    int colFrameSize;
    int extR;
    int extC;
    int rSize;
    int cSize;
    int** zeroFramedAry;
    int** morphAry;
    int** tempAry;
    int** structAry;

    morphology(ifstream& in, ifstream& struc) {
        in >> numImgR >> numImgC >> imgMin >> imgMax;
        struc >> numStructR >> numStructC >> structMin >> structMax >>
            rowOri >> colOri;

        rowFrameSize = numStructR / 2;
        colFrameSize = numStructC / 2;
        extR = rowFrameSize * 2;
        extC = colFrameSize * 2;
        rSize = numImgR + extR;
        cSize = numImgC + extC;

        zeroFramedAry = new int* [rSize];
        morphAry = new int* [rSize];
        tempAry = new int* [rSize];
        structAry = new int* [numStructR];
        for (int i = 0; i < rSize; i++) {
            zeroFramedAry[i] = new int[cSize];
            morphAry[i] = new int[cSize];
            tempAry[i] = new int[cSize];
        }
        for (int i = 0; i < numStructR; i++) {
            structAry[i] = new int[numStructC];
        }
    }

    void zero2DAry(int** ary, int r, int c) {
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                ary[i][j] = 0;
            }
        }
    }

    void loadImg(ifstream& in) {
        for (int i = rowOri; i < rSize - rowFrameSize; i++) {
            for (int j = colOri; j < cSize - colFrameSize; j++) {
                in >> zeroFramedAry[i][j];
            }
        }
    }

    void loadStruct(ifstream& s) {
        for (int i = 0; i < numStructR; i++) {
            for (int j = 0; j < numStructC; j++) {
                s >> structAry[i][j];
            }
        }
    }

    void computeDilation(int** inAry, int** outAry) {
        int i = rowFrameSize;
        while (i < rSize) {
            int j = colFrameSize;
            while (j < cSize) {
                if (inAry[i][j] > 0) {
                    onePixelDilation(i, j, inAry, outAry);
                }
                j++;
            }
            i++;
        }
    }
    void computeErosion(int** inAry, int** outAry) {
        int i = rowFrameSize;
        while (i < rSize) {
            int j = colFrameSize;
            while (j < cSize) {
                if (inAry[i][j] > 0) {
                    onePixelErosion(i, j, inAry, outAry);
                }
                j++;
            }
            i++;
        }
    }

    void onePixelDilation(int i, int j, int** inAry, int** outAry) {
        int iOffset = i - rowOri;
        int jOffset = j - colOri;
        int rIndex = 0;
        while (rIndex < numStructR) {
            int cIndex = 0;

            while (cIndex < numStructC) {
                if (structAry[rIndex][cIndex] > 0) {
                    outAry[iOffset + rIndex][jOffset + cIndex] = 1;
                }
                cIndex++;
            }
            rIndex++;
        }
    }

    void onePixelErosion(int i, int j, int** inAry, int** outAry) {
        int iOffset = i - rowOri;
        int jOffset = j - colOri;
        bool matchFlag = true;
        int rIndex = 0;
        while (matchFlag && rIndex < numStructR) {
            int cIndex = 0;

            while (matchFlag && cIndex < numStructC) {
                
                if (structAry[rIndex][cIndex] > 0 &&
                    inAry[iOffset + rIndex][jOffset + cIndex] <= 0) {
                    matchFlag = false;
                }
                cIndex++;
            }
            rIndex++;
        }
        if (matchFlag) {
            outAry[i][j] = 1;
        }
        else {
            outAry[i][j] = 0;
        }
    }

    void computeClosing(int** zAry, int** mAry, int** tAry) {
        computeDilation(zAry, tAry);
        computeErosion(tAry, mAry);
    }

    void computeOpening(int** zAry, int** mAry, int** tAry) {
        computeErosion(zAry, tAry);
        computeDilation(tAry, mAry);
    }

    void AryToFile(int** ary, ofstream& out) {
        out << numImgR << " " << numImgC << " " << imgMin << " " << imgMax << endl;
        for (int i = rowFrameSize; i < rSize - rowFrameSize; i++) {
            for (int j = colFrameSize; j < cSize - colFrameSize; j++) {
                out << ary[i][j]<<" ";
            }
            out << endl;
        }
    }

    void prettyPrint(int** ary, ofstream& out) {
        out << numImgR << " " << numImgC << " " << imgMin << " " << imgMax << endl;
        for (int i = rowFrameSize ; i < rSize - rowFrameSize; i++) {
            for (int j = colFrameSize; j < cSize - colFrameSize; j++) {
                if (ary[i][j] == 0) {
                    out << ". ";
                }
                else {
                    out << "1 ";
                }
            }
            out << endl;
        }
    }

    void prettyPrintStru(int** ary, ofstream& out) {
        out << numStructR << " " << numStructC << " " << structMin << " " << structMax << endl;
        for (int i = 0; i < numStructR; i++) {
            for (int j = 0; j < numStructC; j++) {
                if (ary[i][j] == 0) {
                    out << ". ";
                }
                else {
                    out << "1 ";
                }
            }
            out << endl;
        }
    }

    void objectExtraction(int** zAry, int** mAry, int** tAry,ofstream& out) {
        out << "entering objectExtraction" << endl;
        zero2DAry(mAry, rSize, cSize);
        computeOpening(zAry, mAry, tAry);
        out << "Printing morphAry After Opening" << endl;
        prettyPrint(mAry, out);
    }

    void fillHoles(int** zAry, int** mAry, int** tAry, ofstream& out) {
        out << "entering fillHoles" << endl;
        zero2DAry(mAry, rSize, cSize);
        computeClosing(zAry, mAry, tAry);
        out << "Printing morphAry After closing" << endl;
        prettyPrint(mAry, out);
    }
};
int main(int argc, char* argv[]) {
    ifstream in(argv[1]);
    ifstream struc(argv[2]);
    ifstream struc2(argv[3]);
    ifstream struc3(argv[4]);
    ifstream img1(argv[5]);
    ofstream img2(argv[6]);
    ofstream out(argv[7]);
   ofstream out2(argv[8]);
    morphology m = morphology(in, struc);
 
    m.zero2DAry(m.zeroFramedAry, m.rSize, m.cSize);
    m.loadImg(in);
    out << "Printing ZeroFramedAry" << endl;
    m.prettyPrint(m.zeroFramedAry, out);

    m.zero2DAry(m.structAry, m.numStructR, m.numStructC);
    m.loadStruct(struc);
    out << "Printing structAry" << endl;
    m.prettyPrintStru(m.structAry, out);

    m.zero2DAry(m.morphAry, m.rSize, m.cSize);
    m.computeDilation(m.zeroFramedAry, m.morphAry);
    out << "Printing morphAry After Dilation" << endl;
    m.prettyPrint(m.morphAry, out);

    m.zero2DAry(m.morphAry, m.rSize, m.cSize);
    m.computeErosion(m.zeroFramedAry, m.morphAry);
    out << "Printing morphAry After Erosion" << endl;
    m.prettyPrint(m.morphAry, out);

    m.zero2DAry(m.morphAry, m.rSize, m.cSize);
    m.computeOpening(m.zeroFramedAry, m.morphAry, m.tempAry);
    out << "Printing morphAry After Opening" << endl;
    m.prettyPrint(m.morphAry, out);

    m.zero2DAry(m.morphAry, m.rSize, m.cSize);
    m.computeClosing(m.zeroFramedAry, m.morphAry, m.tempAry);
    out << "Printing morphAry After Closing" << endl;
    m.prettyPrint(m.morphAry, out);

    morphology t = morphology(img1, struc2);

    out2 << "entering Task1- Object Extraction" << endl;
    t.zero2DAry(t.zeroFramedAry, t.rSize, t.cSize);
    t.loadImg(img1);
    out2 << "Printing ZeroFramedAry" << endl;
    t.prettyPrint(t.zeroFramedAry, out2);

    t.zero2DAry(t.structAry, t.numStructR, t.numStructC);
    t.loadStruct(struc2);
    out2 << "Printing structAry" << endl;
    t.prettyPrintStru(t.structAry, out2);

    t.objectExtraction(t.zeroFramedAry, t.morphAry, t.tempAry, out2);

    t.AryToFile(t.morphAry, img2);
    img2.close();
    ifstream img3(argv[6]);

    morphology t1 = morphology(img3, struc3);

    out2 << "entering Task1- Fill Holes" << endl;
    t1.zero2DAry(t1.zeroFramedAry, t1.rSize, t1.cSize);
    t1.loadImg(img3);
    out2 << "Printing ZeroFramedAry" << endl;
    t1.prettyPrint(t1.zeroFramedAry, out2);

    t1.zero2DAry(t1.structAry, t1.numStructR, t1.numStructC);
    t1.loadStruct(struc3);
    out2 << "Printing structAry" << endl;
    t1.prettyPrintStru(t1.structAry, out2);

    t1.fillHoles(t1.zeroFramedAry, t1.morphAry, t1.tempAry, out2);

    out.close();
    out2.close();
}