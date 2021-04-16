/** ************************************************************************

                Netpbm file handling demo

        uses functions adapted from netpbm-10.73.35
         http://sourceforge.net/projects/netpbm

                    R. Marshall
                     April 2021

compile: g++ -o read --std=c++17 read.cpp
run: ./read [inputfile]

************************************************************************* **/
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdint>
#include <climits>
#include <cerrno>
#include <cstring>
#include <iostream>
#include <vector>
#include <array>
#include <tuple>
#include <fstream>
#include <variant>
#include <omp.h>
#include <stdio.h>

#include "utils.h"

#define MAGIC 0x50 // ascii 'P'
//#define PGMA 0x32  // ascii '2'
//#define PGMB 0x35  // ascii '5'
//#define PBMA 0x31  // ascii '1'
//#define PBMB 0x34  // ascii '4'
//#define PPMA 0x33  // ascii '3'
//#define PPMB 0x34  // ascii '6'

#define PBM_MAGIC1 'P'
#define PBM_MAGIC2 '1'
#define RPBM_MAGIC2 '4'
#define PBM_FORMAT (PBM_MAGIC1 * 256 + PBM_MAGIC2)
#define RPBM_FORMAT (PBM_MAGIC1 * 256 + RPBM_MAGIC2)
#define PBM_TYPE PBM_FORMAT
#define PGM_OVERALLMAXVAL 65535
#define PGM_MAXMAXVAL 255

//#define DEFAULT_MAXVAL 255
#define DEFAULT_INPUT "./input/images/x1b.pgm"

typedef unsigned int gray;
typedef std::vector<std::vector<gray>> pixeldata;
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
/*
  from lib/util/mallocvar.h
*/
static inline void
mallocProduct(unsigned char **      const resultP, 
              unsigned int const factor1,
              unsigned int const factor2) {
/*----------------------------------------------------------------------------
   malloc a space whose size in bytes is the product of 'factor1' and
   'factor2'.  But if that size cannot be represented as an unsigned int,
   return NULL without allocating anything.  Also return NULL if the malloc
   fails.

   If either factor is zero, malloc a single byte.

   Note that malloc() actually takes a size_t size argument, so the
   proper test would be whether the size can be represented by size_t,
   not unsigned int.  But there is no reliable indication available to
   us, like UINT_MAX, of what the limitations of size_t are.  We
   assume size_t is at least as expressive as unsigned int and that
   nobody really needs to allocate more than 4GB of memory.
-----------------------------------------------------------------------------*/
    if (factor1 == 0 || factor2 == 0)
        *resultP = (unsigned char*)malloc(1);
    else {
        if (UINT_MAX / factor2 < factor1) 
            *resultP = NULL;
        else 
            *resultP = (unsigned char*)malloc(factor1 * factor2); 
    }
}

/*
  from lib/util/mallocvar.h
*/
#define MALLOCARRAY(arrayName, nElements) do { \
    unsigned char * array; \
    mallocProduct(&array, nElements, sizeof(arrayName[0])); \
    arrayName = array; \
} while (0)


struct ppm_pixel {
  std::array<int,3> pixel_array;
  static std::string to_string(const ppm_pixel& _ppm_pixel);
};

inline std::string ppm_pixel::to_string(const ppm_pixel& _ppm_pixel) {
  return std::string(std::to_string(_ppm_pixel.pixel_array[0]) + std::string(" ") + std::to_string(_ppm_pixel.pixel_array[1]) + std::string(" ") + std::to_string(_ppm_pixel.pixel_array[2]));
}


struct header_t {
  int size;
  char format[3] = {MAGIC,'\0','\0'};
};

struct pgm_header_t : header_t {

  int width;
  int height;
  gray maxval;

};
//--------------------------------------------------------------------------
// prototypes
//--------------------------------------------------------------------------


// from lib/libpgm1.c
static void
readRpgmRow(FILE * const fileP, gray * const grayrow, pgm_header_t *header) ;


// from lib/pmfileio.c
int get_magic_bytes(FILE * const ifP) ;


// from lib/fileio.c
char pm_getc(FILE * const fileP) ;


// from lib/pmfileio.c
unsigned int m_getuint(FILE * const ifP) ;

// from lib/libpgm1.c
static void
validateComputableSize(unsigned int const cols,
                       unsigned int const rows) ;



/*
  from lib/libpgm1.c
*/
void
pgm_readpgminitrest(FILE * const fileP, 
                    pgm_header_t *header) ;

/*
  from lib/libpgm1.c
*/
static void
readRpgmRow(FILE * const fileP,
            gray * const grayrow, 
            pgm_header_t *header) ;



struct pgm_header_t get_header(FILE * const ifP, char magic) ;


// TODO: Define
// My prototypes!
// All convolution functions have side effects on the vector image. The
// vector is passed by reference and operated on directly to save time from 
// memory allocation (specifically in horizontal_flip and vertical_flip).

bool write_image_to_file(const char format, const std::vector<std::vector<gray>>& image, const std::string& path);

bool write_image_to_file(const char format, const std::vector<std::vector<ppm_pixel>>& image, const std::string& path);

std::tuple<int,int> get_width_height(const std::string& line);


void pm_error(std::string errstr) {
  std::cerr << errstr << std::endl;
}

unsigned char
pm_getrawbyte(FILE * const file) {
    int iby;

    iby = getc(file);
    if (iby == EOF)
        pm_error("EOF / read error reading a one-byte sample");
    return (unsigned char) iby;
}

typedef unsigned char bit;
static bit getbit (FILE * const file) {
    char ch;

    do {
        ch = pm_getc( file );
    } while ( ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' );

    if ( ch != '0' && ch != '1' )
        pm_error( "junk in file where bits should be" );
    
    return ( ch == '1' ) ? 1 : 0;
}

void pbm_readpbmrow( FILE * const file,
                bit * const bitrow,
                int const cols,
                int const format) {

    int col, bitshift;

    switch ( format )
    {
    case PBM_FORMAT:
    for ( col = 0; col < cols; ++col )
        bitrow[col] = getbit( file );
    break;

    case RPBM_FORMAT: {
        unsigned char item;
        bitshift = -1;  item = 0;  /* item's value is meaningless here */
        for ( col = 0; col < cols; ++col ) {
              if ( bitshift == -1 ) {
                    item = pm_getrawbyte( file );
                    bitshift = 7;
                }
              bitrow[col] = ( item >> bitshift ) & 1;
              --bitshift;
          }
    }
    break;

    default:
    pm_error( "can't happen" );
    }
}


//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

inline void print_pixeldata(struct pgm_header_t hh, pixeldata data) {
    for (int r = 0; r < hh.height; ++r){
      for (int c = 0; c < hh.width; ++c){
        std::cerr << data[r][c] << " " ;
      }
      std::cerr << std::endl;
    }
}

/*
  new implementation here, because so many calls to pm_error in 
    the functions from netpbm.  There is another variant that takes 
    more arguments and behaves similarly to fprintf or sprintf.

*/
//------------------------------------------------------------------------

int get_magic_bytes(FILE * const ifP) {

    int ich1, ich2;

    ich1 = getc(ifP);
    ich2 = getc(ifP);
    if (ich1 == EOF || ich2 == EOF)
        std::cerr << "Error reading magic number from Netpbm image stream.  " <<
                  "Most often, this " <<
                  "means your input file is empty." << std::endl;

    return ich1 * 256 + ich2;
} //------------------------------------------------------------------------

char pm_getc(FILE * const fileP) {
  int ich;
  char ch;

  ich = getc(fileP);
  if (ich == EOF)
      pm_error("EOF / read error reading a byte");
  ch = (char) ich;
  
  if (ch == '#') {
      do {
          ich = getc(fileP);
          if (ich == EOF)
              pm_error("EOF / read error reading a byte");
          ch = (char) ich;
      } while (ch != '\n' && ch != '\r');
  }
  return ch;
} //------------------------------------------------------------------------

unsigned int pm_getuint(FILE * const ifP) {
/*----------------------------------------------------------------------------
   Read an unsigned integer in ASCII decimal from the file stream
   represented by 'ifP' and return its value.

   If there is nothing at the current position in the file stream that
   can be interpreted as an unsigned integer, issue an error message
   to stderr and abort the program.

   If the number at the current position in the file stream is too
   great to be represented by an 'int' (Yes, I said 'int', not
   'unsigned int'), issue an error message to stderr and abort the
   program.
-----------------------------------------------------------------------------*/
  char ch;
  unsigned int i;

  do {
      ch = pm_getc(ifP);
  } while (ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r');

  if (ch < '0' || ch > '9')
      pm_error("junk in file where an unsigned integer should be");

  i = 0;
  do {
      unsigned int const digitVal = ch - '0';

      if (i > INT_MAX/10)
          pm_error("ASCII decimal integer in file is "
                   "too large to be processed.  ");
      
      i *= 10;

      if (i > INT_MAX - digitVal)
          pm_error("ASCII decimal integer in file is "
                   "too large to be processed.  ");

      i += digitVal;

      ch = pm_getc(ifP);
  } while (ch >= '0' && ch <= '9');

  return i;
} //------------------------------------------------------------------------

static void
validateComputableSize(unsigned int const cols,
                       unsigned int const rows) {
/*--------------------------------------------------------------------------
   Validate that the dimensions of the image are such that it can be
   processed in typical ways on this machine without worrying about
   overflows.  Note that in C, arithmetic is always modulus
   arithmetic, so if your values are too big, the result is not what
   you expect.  That failed expectation can be disastrous if you use
   it to allocate memory.

   It is very normal to allocate space for a pixel row, so we make sure
   the size of a pixel row, in bytes, can be represented by an 'int'.

   A common operation is adding 1 or 2 to the highest row or
   column number in the image, so we make sure that's possible.
---------------------------------------------------------------------------*/
    if (cols > INT_MAX / (sizeof(gray)) || cols > INT_MAX - 2)
        fprintf(stderr, "image width (%u) too large to be processed", cols);
    if (rows > INT_MAX - 2)
        fprintf(stderr, "image height (%u) too large to be processed", rows);
} //------------------------------------------------------------------------

void
pgm_readpgminitrest(FILE * const fileP, 
                    pgm_header_t *header) {
  gray maxval;
  
  /* Read size. */
  header->width = (int)pm_getuint(fileP);
  header->height = (int)pm_getuint(fileP);

  /* Read maxval. */
  maxval = pm_getuint(fileP);
  if (maxval > PGM_OVERALLMAXVAL)
      fprintf(stderr,
        "maxval of input image (%u) is too large. \
        The maximum allowed by PGM is %u.", 
               maxval, PGM_OVERALLMAXVAL);
  if (maxval == 0)
      pm_error("maxval of input image is zero.");

  header->maxval = maxval;
} //------------------------------------------------------------------------

struct pgm_header_t get_header(FILE * const ifP, char magic) {
  struct pgm_header_t header;

  switch(magic) {
    case '2':
    case '5': {
      std::cerr << "-- read pgm header --" << std::endl;
      pgm_readpgminitrest(
        ifP, &header) ;
      
      break;
    }
    default: {
      //std::cerr << "add more cases to handle more types" << std::endl;
    }
  }

  validateComputableSize(header.width, header.height);

  return header;
} //------------------------------------------------------------------------

static void
readRpgmRow(FILE * const fileP,
            gray * const grayrow, 
            pgm_header_t *header) {

    unsigned int const bytesPerSample = header->maxval < 256 ? 1 : 2;
    int          const bytesPerRow    = header->width * bytesPerSample;
    
    unsigned char * rowBuffer;
    char * error;// = nullptr;

    MALLOCARRAY(rowBuffer, bytesPerRow);
    if (rowBuffer == NULL)
        sprintf(error, "Unable to allocate memory for row buffer "
                    "for %u columns", header->width);
    else {
      size_t rc;
      rc = fread(rowBuffer, 1, bytesPerRow, fileP);
      if (rc == 0)
        sprintf(error, "Error reading row.  fread() errno=%d (%s)",
                      errno, strerror(errno));
      else if (rc != bytesPerRow)
        sprintf(error, "Error reading row.  Short read of %u bytes "
                      "instead of %u", (unsigned)rc, bytesPerRow);
      else {
        if (header->maxval < 256) {
            unsigned int col;
            for (col = 0; col < header->width; ++col){
                grayrow[col] = (gray)rowBuffer[col];
              }

        } else { // we won't deal with wide formats    
        }
        if (header->maxval == 255 || header->maxval == 65535) {
            /* There's no way a sample can be invalid, so we don't need to look at
               the samples individually.
            */

        } else {
          unsigned int col;
          for (col = 0; col < header->width; ++col) {
            if (grayrow[col] > header->maxval) {
                 fprintf(stderr, "gray value %u is greater than maxval (%u)", grayrow[col], header->maxval );

                return;
            }
          }
        }
      }
      free(rowBuffer);
    }
    if (error) {
      std::cerr << "error in binary file read" << std::endl;
    }
} //------------------------------------------------------------------------




// ====================== MAIN =============================





int main (int argc, char *argv[]) {
  double start(omp_get_wtime());
  std::string filename = DEFAULT_INPUT;

  if (argc > 1) {
    filename = std::string(argv[1]);
  }

  pgm_header_t header;
  FILE* imagefile;
  imagefile = fopen(filename.c_str(), "r"); // open for reading, ascii mode for header

  /**
     magic number really only uses the lowest 16 bits.  Example: x1b.pgm magic number 
     looks like this in binary:

     00000000000000000101000000110101

     split into bytes, it looks like:

     00000000 00000000 01010000 00110101

     cast to a char, it will only retain the last 8 bits to get the '1' thru '6'
     right shift 8 bits and cast to a char (optionally here).  It should be 01010000 ('P')
     for each case, you can use this part to verify

  */
  int k = get_magic_bytes(imagefile); 
  char mn = (char)k; // cast the lowest 8 bits of k into a char
  pgm_header_t hh = get_header(imagefile, mn);
  hh.format[1] = (char)mn;
  /*
  std::cerr << "image file:    " << filename << std::endl;
  std::cerr << "magic byte hi: " << (char)(k>>8) << std::endl;
  std::cerr << "magic byte lo: " << mn << std::endl;
  std::cerr << "format:        " << hh.format << std::endl;
  std::cerr << "width:         " << hh.width << std::endl;
  std::cerr << "height:        " << hh.height << std::endl;
  std::cerr << "maxval:        " << hh.maxval << std::endl;
  */ 
  gray *grayrow;

  // std::variant of two types of vectors.. depending on which case is used in
  // the switch statement, it could be std::vector<std::vector<gray>> or 
  // std::vector<std::vector<ppm_pixel>>. 
  // Instead of defining two objects for each, you can store both in the same
  // portion of memory with std::variant (much like a union). 
  // This is useful when checking the types with std::holds_alternative and 
  // calling templated functions (provided in utils.h) to perform the
  // convolution operations on images of different types.
  std::variant<std::vector<std::vector<gray>>,std::vector<std::vector<ppm_pixel>>> image;
  int max_value{}; 
  bool closed = false;

  switch (hh.format[1]) {

    case '1': {
      // Case 1: format of pbm ascii
      std::cerr << "format:        pbm ascii" << std::endl;
      // Closing the imagefile as we will use a std::ifstream instead
      fclose(imagefile);
      closed = true;
      std::ifstream input_file;
      input_file.open(filename);
      std::string line;
      int width, height;
      std::getline(input_file,line);
      if (!(line == "P1")) {
        std::cerr << "Not a pbm ascii file" << std::endl;
        return -1;
      }
      while (true) {
        std::getline(input_file,line);
        // Ignoring # comments
        if (line.find("#") != std::string::npos)
          continue;
        std::tie(width, height) = get_width_height(line);
        // Ensuring we get proper values for width and height
        // If not, return -1
        if (width == 0 && height == 0) {
          return -1;
        }
        break;
      }
      // Resizing 2d vector to height x width (rows x columns)
      std::vector<std::vector<gray>> temp_image(height, std::vector<gray>(width));
      // Reading in contents of input_file into 2d vector
      // White spaces are the delimiters for >>, so we do not 
      // need to worry about spaces or newline characters.
      for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
          input_file >> temp_image[i][j]; 
        }
      }

      // moving temp_image into std::variant image
      image = std::move(temp_image);
      break;  }


    case '2': {
      std::vector<std::vector<gray>> temp_image;
      std::cerr << "format:        pgm ascii" << std::endl;

      grayrow = (gray*)malloc(hh.width * sizeof(gray));

      for (int row = 0; row < hh.height; ++row) {
        for (gray col = 0; col < (unsigned)hh.width; ++col) {
          grayrow[col] = pm_getuint(imagefile);

          if (grayrow[col] > hh.maxval) {
              fprintf(stderr, "value out of bounds (%u > %u)",
                       grayrow[col], hh.maxval);
          }
        }
        std::vector<gray> pgmrow(grayrow, grayrow+hh.width);
        temp_image.push_back(pgmrow);
      }
      // uncomment the line below to print the data section
      // print_pixeldata(hh, image);
      free(grayrow);
      // moving temp_image into std::variant image
      image = std::move(temp_image);

      break;  }


    case '3': {
      std::cerr << "format:        ppm ascii" << std::endl;
      // Closing the imagefile as we will use a std::ifstream instead
      fclose(imagefile);
      closed = true;
      std::ifstream input_file;
      input_file.open(filename);
      std::string line;
      int width, height;
      std::getline(input_file,line);
      if (!(line == "P3")) {
        std::cerr << "Not a ppm ascii file" << std::endl;
        return -1;
      }
      while (true) {
        std::getline(input_file,line);
        // Ignoring # comments
        if (line.find("#") == !std::string::npos)
          continue;
        std::tie(width, height) = get_width_height(line);
        // Ensuring we get proper values for width and height
        // If not, return -1
        if (width == 0 && height == 0) {
          return -1;
        }
        std::getline(input_file,line);
        break;
      }
      // Resizing 2d vector to height x width (rows x columns)
      std::vector<std::vector<ppm_pixel>> temp_image(height, std::vector<ppm_pixel>(width));
      // Reading in contents of input_file into 2d vector
      // White spaces are the delimiters for >>, so we do not 
      // need to worry about spaces or newline characters.
      for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
          input_file >> temp_image[i][j].pixel_array[0];
          input_file >> temp_image[i][j].pixel_array[1];
          input_file >> temp_image[i][j].pixel_array[2];
        }
      }
      // moving temp_image into std::variant image
      image = std::move(temp_image);
      break;  }

    /*
    case '4': {
      
      std::cerr << "format:        pbm binary" << std::endl;
      // Case 4: format of pbm binary 
      // Closing the imagefile as we will use a std::ifstream instead
       
      fclose(imagefile);
      std::ifstream input_file;
      input_file.open(filename);
      std::string line;
      std::getline(input_file,line);
      std::stringstream ss;
      std::cout << line << std::endl;
      ss << line;
      std::array<std::string,3> string_arr;
      for (int i = 0; i < 3; ++i) {
        std::string val;
        ss >> val;
        string_arr[i] = val;
      }
      ss >> string_arr[0];
      ss >> string_arr[1];
      ss >> string_arr[2];
      std::cout << string_arr[0] << std::endl;
      std::cout << string_arr[1] << std::endl;
      std::cout << string_arr[2] << std::endl;
      if (string_arr[0] != "P4")
        return -1;
      int width = std::stoi(string_arr[1]);
      int height = std::stoi(string_arr[2]);
      if (width == 0 && height == 0) 
        return -1;
      // Switch to read binary
      input_file.close();
      /*
      input_file.open(filename, std::ios::binary);
      std::getline(input_file,line);
      std::cout << "I am here" << std::endl;
      // Resizing 2d vector to height x width (rows x columns)
      std::vector<std::vector<gray>> temp_image(height, std::vector<gray>(width));
      // Reading in contents of input_file into 2d vector
      // White spaces are the delimiters for >>, so we do not 
      // need to worry about spaces or newline characters.
      std::bitset<1> v;
      std::bitset<1>::reference ref = v[0];
      for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
          input_file.read(reinterpret_cast<char*>(&ref), 1);
          std::cout << v.to_string() << ' ';
          temp_image[i][j] = v.to_ulong();
        }
        std::cout << '\n';
      }
      

      return -1;
      // moving temp_image into std::variant image
      
      image = std::move(temp_image);
      
      std::vector<std::vector<gray>> temp_image;
      grayrow = (gray*)
      grayrow = (gray*)malloc(width * sizeof(gray));
      for (int row = 0; row < height; ++row) {
        readPbmRow(imagefile, grayrow, &hh);
        pbm_readpbmrow(imagefile);
        std::vector<gray> pgmrow(grayrow, grayrow+hh.width);
        temp_image.push_back(pgmrow);
      }

      break;  }





    case '5': {
      std::cerr << "format:        pgm binary" << std::endl;
      std::vector<std::vector<gray>> temp_image;
      grayrow = (gray*)malloc(hh.width * sizeof(gray));

      for (int row = 0; row < hh.height; ++row) {
        readRpgmRow(imagefile, grayrow, &hh);
        std::vector<gray> pgmrow(grayrow, grayrow+hh.width);
        temp_image.push_back(pgmrow);
      }
      // uncomment the line below to print the data section
      // print_pixeldata(hh, image);
      free(grayrow);
      image = std::move(temp_image);
      
      break;  }


    case '6': {
      std::cerr << "format:    ppm binary" << std::endl;
      break;  }

    */ 
    default: {
      std::cerr << "unknown file type/format" << std::endl;
      return -1;
    }
  }



  // ================== 
  // Image Convolution 
  // At this point we should have our png file filled inside the image
  // vector;
  // ================== 
  

  double parallel_time;
  unsigned short n_threads = static_cast<unsigned short>(atoi(argv[3]));

  switch (atoi(argv[2])) {
    // 90 degrees clockwise
    case 0 :
      if (std::holds_alternative<std::vector<std::vector<gray>>>(image))
        parallel_time = ninety_degree_flip(std::get<0>(image), n_threads);
      else if (std::holds_alternative<std::vector<std::vector<ppm_pixel>>>(image))
        parallel_time = ninety_degree_flip(std::get<1>(image), n_threads);
      break;
    
    // vertical flip
    case 1 : 
      if (std::holds_alternative<std::vector<std::vector<gray>>>(image))
        parallel_time = vertical_flip(std::get<0>(image), n_threads);
      else if (std::holds_alternative<std::vector<std::vector<ppm_pixel>>>(image))
        parallel_time = vertical_flip(std::get<1>(image), n_threads);
      break;

    // horizontal flip
    case 2 :
       if (std::holds_alternative<std::vector<std::vector<gray>>>(image))
        parallel_time = horizontal_flip(std::get<0>(image), n_threads);
      else if (std::holds_alternative<std::vector<std::vector<ppm_pixel>>>(image))
        parallel_time = horizontal_flip(std::get<1>(image), n_threads);
      break;
    
    default : 
      break;
  }
    
  std::stringstream ss;
  if (std::holds_alternative<std::vector<std::vector<gray>>>(image)) { 
    ss << filename.substr(filename.find_last_of('/')+1, filename.find_last_of('.') - filename.find_last_of('/') - 1)  << '_' << atoi(argv[2]);
    if(!write_image_to_file(hh.format[1], std::get<0>(image), ss.str())) { 
      std::cerr << "Error in writing to file" << std::endl;
      return -1;
    }

  } 
  else if (std::holds_alternative<std::vector<std::vector<ppm_pixel>>>(image)) {
    ss << filename.substr(filename.find_last_of('/')+1,  filename.find_last_of('.') - filename.find_last_of('/') - 1) << '_' << atoi(argv[2]);
    if(!write_image_to_file(hh.format[1], std::get<1>(image), ss.str())) {
      std::cerr << "Error in writing to file" << std::endl;
      return -1;
    }
  }
  
  // TODO Need a conditional to set if this is acually close 
  if (!closed) {
    fclose(imagefile);
  }

  double end(omp_get_wtime());
  double total_time(end - start);
  

  // Printing info
  if (std::holds_alternative<std::vector<std::vector<gray>>>(image)) {
    print_info(argv, std::get<0>(image), n_threads, parallel_time, total_time);
  }
  else if (std::holds_alternative<std::vector<std::vector<ppm_pixel>>>(image)) {
    print_info(argv, std::get<1>(image), n_threads, parallel_time, total_time); 
  }
  
  return 0;
}


// Write image to file function for std::vector<std::vector<gray>>&;
// corresponds to pbm or pgm file types. 
bool write_image_to_file(const char format, const std::vector<std::vector<gray>>& image, const std::string& path) {
  switch (format) {
    case '1': {
      // pbm ascii
      std::ofstream output_file;
      output_file.open(path+std::string(".pbm"));
      if (!output_file.is_open())
        return false;
      output_file << "P1\n" << image[0].size() << ' ' << image.size() << '\n';
      for (int i = 0; i < image.size(); ++i) {
        for (int j = 0; j < image[0].size(); ++j) {
          output_file << image[i][j] << ' ';
        }
        output_file << '\n';
      }
      output_file.close();

      break;}

    case '2': {
      // pgm ascii
      std::ofstream output_file;
      output_file.open(path+std::string(".pgm"));
      if (!output_file.is_open())
        return false;
      output_file << "P2\n" << image[0].size() << ' ' << image.size() << '\n' << 255 << '\n';
      for (int i = 0; i < image.size(); ++i) {
        for (int j = 0; j < image[0].size(); ++j) {
          output_file << image[i][j] << ' ';
        }
        output_file << '\n';
      }
      output_file.close();
      break;}
    /*
    case '4': {
      // pbm binary
      break;}

    case '5': {
      // pgm binary
      break;}
    */
    default: {
      std::cerr << "Unsupported file format to write to" << std::endl;
      return false;
      break;}
  } 
  return true;
}


// Overloaded write_image_to_file() function that takes in
// std::vector<std::vector<ppm_pixel>>& image.
bool write_image_to_file(const char format, const std::vector<std::vector<ppm_pixel>>& image, const std::string& path) {
  
  // Only possible cases are '3' and '6' since this is an overloaded function
  // that takes in ppm image vector.
  switch(format) {
    // ppm ascii
    case '3' : { 
      std::ofstream output_file;
      output_file.open(path+std::string(".ppm"));
      if (!output_file.is_open())
        return false;
      output_file << "P3\n" << image[0].size() << ' ' << image.size() << '\n' << 255 << '\n';
      std::stringstream ss;
      std::string val;
      for (int i = 0; i < image.size(); ++i) {
        for (int j = 0; j < image[0].size(); ++j) {
          // Inline static member function of stuct ppm_pixel
          // Defined on line 97
          val = ppm_pixel::to_string(image[i][j]);
          // Since the ppm standard does not allow for line sizes larger
          // than 70 characters, I needed to use this if statement to 
          // check if concatenating val to ss would result in a string
          // with a size larger than 70.
          // If so, print contents of ss, clear ss, and write val to ss.
          // Else, write val to ss.
          if ((ss.gcount() + val.size()) > 70) {
            ss << '\n';
            output_file << ss.str();
            ss.str(std::string());
            ss << val << ' ';
          } else {
            ss << val << ' ';
          }
        }
        output_file << ss.str() << '\n';
        ss.str(std::string());
      }
      output_file << ss.str();
      output_file.close();
    break; }
    /*  
    // ppm binary
    case '6' : {

    break;}
    */

    default : {
  
    break;}
  }
  return true;
}


// Parses line and obtains width and height.
// Returns tuple to std::tie into values.
// If numeric values cannot be deduced, std::invalid_argument will be thrown.
// In that case, a tuple of contents {0,0} will be returned to signify error. 
std::tuple<int,int> get_width_height(const std::string& line) {
  try {
    return std::make_tuple(std::stoi(line.substr(0, line.find(' '))), std::stoi(line.substr(line.find(' '))));
  }
  catch(std::invalid_argument) {
    std::cerr << "Could not convert width and height to numeric type" << std::endl;
    return std::make_tuple(0,0);
  }
}


