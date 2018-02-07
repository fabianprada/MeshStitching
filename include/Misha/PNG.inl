#include <stdio.h>
#include <vector>
#ifdef _WIN32
#include <PNG/png.h>
#else
#include <libpng/png.h>
#endif


unsigned char* PNGReadColor( const char* fileName , int& width , int& height )
{
	png_structp png_ptr =
		png_create_read_struct(PNG_LIBPNG_VER_STRING,
		0, // (png_voidp)user_error_ptr
		0, // user_error_fn
		0  // user_warning_fn
		);
	if(!png_ptr)	return NULL;
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if(!info_ptr)	return NULL;
	png_infop end_info = png_create_info_struct(png_ptr);
	if(!end_info)	return NULL;

	FILE* fp = fopen( fileName , "rb" );
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open file for reading: %s\n" , fileName ) , exit( 0 );
	png_init_io( png_ptr , fp );
	unsigned char* pixels;


	// tell the library if we have already read any bytes from header
	// png_set_sig_bytes(png_ptr, 0);
	// callback to handle user chunk data
	// png_set_read_user_chunk_fn(png_ptr,user_chunk_ptr,read_chunk_callback);
	// callback used to control a progress meter
	// png_set_read_status_fn(png_ptr, read_row_callback);
	//
	{
		// high-level read
		int png_transforms=(PNG_TRANSFORM_STRIP_16  | // 16-bit to 8-bit
			PNG_TRANSFORM_PACKING   | // expand 1, 2, and 4-bit
			0);
		png_read_png(png_ptr, info_ptr, png_transforms, NULL);
		width=png_get_image_width(png_ptr,info_ptr);
		height=png_get_image_height(png_ptr,info_ptr);
		int ncomp=png_get_channels(png_ptr,info_ptr);
		int bit_depth=png_get_bit_depth(png_ptr,info_ptr);
		int color_type= png_get_color_type(png_ptr,info_ptr);
		if(width<=0 || height<=0)				return NULL;
		if(ncomp<1 || ncomp>4)					return NULL;
		if(bit_depth!=8)						return NULL;
		if(color_type==PNG_COLOR_TYPE_PALETTE && 0)
		{
			png_set_palette_to_rgb(png_ptr);
			png_set_expand(png_ptr);
		}
		png_bytep* row_pointers;   // [height]
		row_pointers = png_get_rows(png_ptr, info_ptr);

		pixels = new unsigned char[ width * height * 3 ];
		if( !pixels ) fprintf( stderr , "[ERROR] Failed to allocate pixels: %d x %d x 3\n" , width , height ) , exit( 0 );
		for(int y=0;y<height;y++)
		{
			unsigned char* buf=(unsigned char*)(row_pointers[y]);
			for(int x=0;x<width;x++)
				for(int z=0;z<ncomp;z++)
					if(color_type==PNG_COLOR_TYPE_PALETTE)
					{
						png_color clr=info_ptr->palette[*buf++];
						pixels[ (y*width+x)*3+0 ] = clr.red;
						pixels[ (y*width+x)*3+1 ] = clr.green;
						pixels[ (y*width+x)*3+2 ] = clr.blue;
					}
					else pixels[ (y*width+x)*3+z ] = *buf++;
		}
	}
	png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
	fclose( fp );
	return pixels;
}
void PNGWriteColor( const char* fileName , const unsigned char* pixels , int width , int height )
{
	FILE* fp = fopen( fileName , "wb" );
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open file for writing: %s\n" , fileName ) , exit( 0 );
	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(!png_ptr)	return;
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if(!info_ptr)	return;
	png_init_io(png_ptr, fp);
	// turn off compression or set another filter
	// png_set_filter(png_ptr, 0, PNG_FILTER_NONE);
	png_set_IHDR(png_ptr, info_ptr, width , height ,
		8,PNG_COLOR_TYPE_RGB,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	if (0) {                    // high-level write
		std::vector<unsigned char> matrix( width * height * 3 );
		std::vector<png_bytep> row_pointers( height );
		for(int y=0;y<height;y++)
		{
			row_pointers[y]=&matrix[y*width*3];
			unsigned char* buf=&matrix[y*width*3];
			for(int x=0;x<width;x++)
				for(int z=0;z<3;z++)
					*buf++ = pixels[ (y*width+x)*3 + z ];
		}
		png_set_rows(png_ptr, info_ptr, &row_pointers[0]);
		int png_transforms=0;
		png_write_png(png_ptr, info_ptr, png_transforms, NULL);
	} else {                    // low-level write
		png_write_info(png_ptr, info_ptr);
		// png_set_filler(png_ptr, 0, PNG_FILLER_AFTER);
		//  but no way to provide GRAY data with RGBA fill, so pack each row
		std::vector<unsigned char> buffer(width*3);
		for(int y=0;y<height;y++)
		{
			unsigned char* buf=&buffer[0];
			for(int x=0;x<width;x++)
				for(int z=0;z<3;z++)
					*buf++ = pixels[ (y*width+x)*3+z ];
			png_bytep row_pointer=&buffer[0];
			png_write_row(png_ptr, row_pointer);
		}
	}
	png_write_end(png_ptr, NULL);
	png_destroy_write_struct(&png_ptr, &info_ptr);
	fclose( fp );
}
