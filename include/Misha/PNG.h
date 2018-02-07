#ifndef PNG_INCLUDED
#define PNG_INCLUDED
void PNGWriteColor( const char* fileName , const unsigned char* pixels , int width , int height );
unsigned char* PNGReadColor( const char* fileName , int& width , int& height );

#include "PNG.inl"
#endif //PNG_INCLUDED
