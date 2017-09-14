//
// Created by nicolet on 14/09/17.
//

#include "Image.h"

//  Image.cpp
//

#ifndef _Image_h
#define _Image_h
#define WIDTH 128
#define HEIGHT 128

#include "Image.h"
#include <cmath>

#endif

using namespace std;
typedef unsigned char unchar;

//Constructor

MImage::MImage(const char* fileName){

    imageData = new unchar* [HEIGHT]; // create new array size: height of image.
    filteredData = new unchar* [HEIGHT];// create new array size: height of image.

    for (int i = 0; i < HEIGHT; i++) {

        imageData[i] = new unchar [WIDTH]; //create matrix.
        filteredData[i] = new unchar [WIDTH]; //create matrix.
    }

    //image I/O
    pInFile = new ifstream;
    pInFile->open(fileName, ios::in | ios::binary); // open fileName and read as binary.
    pInFile->seekg(0, ios::beg); //pos filter at beginning of image file.
    pInFile->read(reinterpret_cast<char*>(imageHeaderData),1078); //read bmp header data into array.

    for(int i=0; i<HEIGHT; i++) {
        pInFile->read(reinterpret_cast<char*>(imageData[i]),WIDTH);//read row into each array entry.
    }

    pInFile->close(); //close stream.

    char m_smoothFilter[3][3] = {
            {1,1,1},
            {1,2,1},
            {1,1,1}
    };

}

unchar** MImage::getData() const {
    return imageData;
}

MImage::~MImage(){

    delete pInFile;
    delete pOutFile;

    for(int i=0; i<HEIGHT; i++){
        delete[] imageData[i];
        delete[] filteredData[i];
    }

    delete[] imageData;
    delete[] filteredData;
}

void MImage::write(const char* fileName) {

    pOutFile = new ofstream;
    pOutFile->open(fileName, ios::out | ios::trunc | ios::binary);
    pOutFile->write(reinterpret_cast<char*>(imageHeaderData), 1078); //write header data onto output

    for(int i = 0; i < HEIGHT; i++){

        pOutFile->write(reinterpret_cast<char*>(filteredData[i]),WIDTH); // write new image data.

    }

    pOutFile->close(); //close stream
}