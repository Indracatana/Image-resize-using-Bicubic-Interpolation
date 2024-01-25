#include <stdio.h>

#include <stdlib.h>

#include <math.h>

#include <pthread.h>

#include "lodepng.cpp"

using namespace std;

struct pthread_struct {
  int id;
  int cores;
  int h, w;
  float factor;
  unsigned char ** imgR;
  unsigned char ** imgG;
  unsigned char ** imgB;
  unsigned char ** img2gR;
  unsigned char ** img2gG;
  unsigned char ** img2gB;
  pthread_mutex_t mutex_R;
  pthread_mutex_t mutex_G;
  pthread_mutex_t mutex_B;
};

// Alocare matrice de dimensiunea h*w
unsigned char ** createMatrix(int h, int w) {
  unsigned char ** mat;
  int i;

  //Alocare si verificarea alocarii la fiecare pas
  mat = (unsigned char ** ) malloc(h * sizeof(unsigned char * ));
  if (mat == NULL) return (NULL);

  mat[0] = (unsigned char * ) malloc(h * w * sizeof(unsigned char));
  if (mat[0] == NULL) return (NULL);

  for (i = 1; i < h; ++i)
    mat[i] = mat[i - 1] + w;
  return (mat);
}

// Functie de minim
int min(int a, int b) {
  return a < b ? a : b;
}

// Functie de calcul a pixelului folosind cei 4x4 pixeli din jur cu bicubic interpolation
double bicubicpol(double x, double y, double p[4][4]) {

  double a00, a01, a02, a03;
  double a10, a11, a12, a13;
  double a20, a21, a22, a23;
  double a30, a31, a32, a33;
  double x2 = x * x;
  double x3 = x2 * x;
  double y2 = y * y;
  double y3 = y2 * y;

  a00 = p[1][1];
  a01 = -.5 * p[1][0] + .5 * p[1][2];
  a02 = p[1][0] - 2.5 * p[1][1] + 2 * p[1][2] - .5 * p[1][3];
  a03 = -.5 * p[1][0] + 1.5 * p[1][1] - 1.5 * p[1][2] + .5 * p[1][3];
  a10 = -.5 * p[0][1] + .5 * p[2][1];
  a11 = .25 * p[0][0] - .25 * p[0][2] - .25 * p[2][0] + .25 * p[2][2];
  a12 = -.5 * p[0][0] + 1.25 * p[0][1] - p[0][2] + .25 * p[0][3] + .5 * p[2][0] - 1.25 * p[2][1] + p[2][2] - .25 * p[2][3];
  a13 = .25 * p[0][0] - .75 * p[0][1] + .75 * p[0][2] - .25 * p[0][3] - .25 * p[2][0] + .75 * p[2][1] - .75 * p[2][2] + .25 * p[2][3];
  a20 = p[0][1] - 2.5 * p[1][1] + 2 * p[2][1] - .5 * p[3][1];
  a21 = -.5 * p[0][0] + .5 * p[0][2] + 1.25 * p[1][0] - 1.25 * p[1][2] - p[2][0] + p[2][2] + .25 * p[3][0] - .25 * p[3][2];
  a22 = p[0][0] - 2.5 * p[0][1] + 2 * p[0][2] - .5 * p[0][3] - 2.5 * p[1][0] + 6.25 * p[1][1] - 5 * p[1][2] + 1.25 * p[1][3] + 2 * p[2][0] - 5 * p[2][1] + 4 * p[2][2] - p[2][3] - .5 * p[3][0] + 1.25 * p[3][1] - p[3][2] + .25 * p[3][3];
  a23 = -.5 * p[0][0] + 1.5 * p[0][1] - 1.5 * p[0][2] + .5 * p[0][3] + 1.25 * p[1][0] - 3.75 * p[1][1] + 3.75 * p[1][2] - 1.25 * p[1][3] - p[2][0] + 3 * p[2][1] - 3 * p[2][2] + p[2][3] + .25 * p[3][0] - .75 * p[3][1] + .75 * p[3][2] - .25 * p[3][3];
  a30 = -.5 * p[0][1] + 1.5 * p[1][1] - 1.5 * p[2][1] + .5 * p[3][1];
  a31 = .25 * p[0][0] - .25 * p[0][2] - .75 * p[1][0] + .75 * p[1][2] + .75 * p[2][0] - .75 * p[2][2] - .25 * p[3][0] + .25 * p[3][2];
  a32 = -.5 * p[0][0] + 1.25 * p[0][1] - p[0][2] + .25 * p[0][3] + 1.5 * p[1][0] - 3.75 * p[1][1] + 3 * p[1][2] - .75 * p[1][3] - 1.5 * p[2][0] + 3.75 * p[2][1] - 3 * p[2][2] + .75 * p[2][3] + .5 * p[3][0] - 1.25 * p[3][1] + p[3][2] - .25 * p[3][3];
  a33 = .25 * p[0][0] - .75 * p[0][1] + .75 * p[0][2] - .25 * p[0][3] - .75 * p[1][0] + 2.25 * p[1][1] - 2.25 * p[1][2] + .75 * p[1][3] + .75 * p[2][0] - 2.25 * p[2][1] + 2.25 * p[2][2] - .75 * p[2][3] - .25 * p[3][0] + .75 * p[3][1] - .75 * p[3][2] + .25 * p[3][3];

  return (a00 + a01 * y + a02 * y2 + a03 * y3) +
    (a10 + a11 * y + a12 * y2 + a13 * y3) * x +
    (a20 + a21 * y + a22 * y2 + a23 * y3) * x2 +
    (a30 + a31 * y + a32 * y2 + a33 * y3) * x3;

}

// Salvare date matrice 1D(vector) in imagine PNG
int savePNG(const char * file, unsigned char * img, int width, int height) {
  unsigned char * png;
  size_t pngsize;
  int error = lodepng_encode24( & png, & pngsize, img, width, height);
  if (!error) {
    lodepng_save_file(png, pngsize, file);
  }

  if (error)
    printf("\tEroare %u la deschiderea imaginii %s: %s\n", error, file,
      lodepng_error_text(error));

  free(png);
  return (error);
}

void * f(void * arg) {
  pthread_struct argument = * (pthread_struct * ) arg;
  int start = argument.id * (double)(argument.factor * argument.h) / argument.cores;
  int end = min((argument.id + 1) * (double)(argument.factor * argument.h) / argument.cores, (argument.factor * argument.h));

  double arrR[4][4] = {
    {
      0
    }
  };
  double arrG[4][4] = {
    {
      0
    }
  };
  double arrB[4][4] = {
    {
      0
    }
  };
  int i, j;

  // Calcul pixeli rezulatati in functie de cei 4x4 pixeli inconjuratori
  for (int i = start; i < end; i++)
    for (int j = 0; j < (int)(argument.factor * argument.w); j++) {
      pthread_mutex_lock( & argument.mutex_R);
      for (int l = 0; l < 4; l++) {
        for (int k = 0; k < 4; k++) {
          arrR[l][k] = (double) argument.imgR[(int)(i / argument.factor)][(int)(j / argument.factor)];
        }
      }
      pthread_mutex_unlock( & argument.mutex_R);

      argument.img2gR[i][j] = (unsigned char) bicubicpol(min(((double) i - ((int)(i / argument.factor) * argument.factor)) / argument.factor, argument.h),
        min(((double) j - ((int)(j / argument.factor) * argument.factor)) / argument.factor, argument.w),
        arrR);
    }
  for (int i = start; i < end; i++)
    for (int j = 0; j < (int)(argument.factor * argument.w); j++) {
      pthread_mutex_lock( & argument.mutex_G);
      for (int l = 0; l < 4; l++) {
        for (int k = 0; k < 4; k++) {
          arrG[l][k] = (double) argument.imgG[(int)(i / argument.factor)][(int)(j / argument.factor)];
        }
      }
      pthread_mutex_unlock( & argument.mutex_G);

      argument.img2gG[i][j] = (unsigned char) bicubicpol(min(((double) i - ((int)(i / argument.factor) * argument.factor)) / argument.factor, argument.h),
        min(((double) j - ((int)(j / argument.factor) * argument.factor)) / argument.factor, argument.w),
        arrG);
    }

  for (int i = start; i < end; i++)
    for (int j = 0; j < (int)(argument.factor * argument.w); j++) {
      pthread_mutex_lock( & argument.mutex_B);
      for (int l = 0; l < 4; l++) {
        for (int k = 0; k < 4; k++) {
          arrB[l][k] = (double) argument.imgB[(int)(i / argument.factor)][(int)(j / argument.factor)];
        }
      }
      pthread_mutex_unlock( & argument.mutex_B);

      argument.img2gB[i][j] = (unsigned char) bicubicpol(min(((double) i - ((int)(i / argument.factor) * argument.factor)) / argument.factor, argument.h),
        min(((double) j - ((int)(j / argument.factor) * argument.factor)) / argument.factor, argument.w),
        arrB);
    }

}

// Citire date din fisier PNG
unsigned char * readPNG(const char * file, unsigned int & width, unsigned int & height,
  unsigned int & bitdepth, unsigned int & bitsXpixel, unsigned int & channels,
  unsigned int & isGrey, unsigned int & haveAlpha) {
  unsigned error;
  unsigned char * image;
  unsigned char * png = 0;
  size_t pngsize;
  LodePNGState state;

  lodepng_state_init( & state);

  error = lodepng_load_file( & png, & pngsize, file);
  if (!error) error = lodepng_decode( & image, & width, & height, & state, png, pngsize);
  if (error) printf("Eroare %u: %s\n", error, lodepng_error_text(error));
  free(png);

  LodePNGColorMode & color = state.info_png.color;

  bitdepth = color.bitdepth;
  bitsXpixel = lodepng_get_bpp( & color);
  channels = lodepng_get_channels( & color);
  isGrey = lodepng_is_greyscale_type( & color);
  haveAlpha = lodepng_can_have_alpha( & color);

  lodepng_state_cleanup( & state);
  return (image);
}

/*  Dintr-o imagine salvata ca matrice 1D, această functie
 *  returneaza o matrice 2D care are componenta de culoare indicata de index:
 *  index = 0  - canal rosu
 *  index = 1  - canal verde
 *  index = 2  - canal albastru
 */
unsigned char ** getChannel2D(unsigned char * img1D, unsigned int h, unsigned int w,
  unsigned int index) {
  unsigned int i, j, k, l;
  unsigned int channels = 3;
  unsigned char ** mat;

  // Alocarea memoriei
  mat = (unsigned char ** ) malloc(h * sizeof(unsigned char * ));
  if (mat == NULL) return (NULL);
  mat[0] = (unsigned char * ) malloc(h * w * sizeof(unsigned char));
  if (mat[0] == NULL) return (NULL);
  for (i = 1; i < h; ++i)
    mat[i] = mat[i - 1] + w;

  // Citirea datelor
  l = (channels + 1);
  for (i = 0; i < h; i++) {
    k = i * w * (channels + 1);
    for (j = 0; j < w; j++) {
      mat[i][j] = img1D[j * l + index + k];
    }
  }

  return (mat);
}

/* Converteste o matrice 2D care contine informatii despre imagine intr-o matrice 1D
 * pentru a putea salva ca imagine PNG
 */
unsigned char * convert2Dto1D(unsigned char ** img2D, unsigned int h, unsigned int w) {
  unsigned char * img1D = (unsigned char * ) malloc(sizeof(unsigned char) * h * w * 3);
  unsigned int i, j, k, l;
  unsigned char val;

  for (i = 0; i < h; i++) {
    k = 3 * w * i;
    for (j = 0; j < w; j++) {
      val = (unsigned char) img2D[i][j];
      l = 3 * j + k;
      img1D[l] = val;
      img1D[l + 1] = val;
      img1D[l + 2] = val;
    }
  }
  return (img1D);
}

// Combina informatiile din trei matrici 2D pentru a forma o matrice 1D pentru o imagine color
unsigned char * convert2Dto1D(unsigned char ** imgR, unsigned char ** imgG,
  unsigned char ** imgB, unsigned int h, unsigned int w) {
  unsigned char * img1D = (unsigned char * ) malloc(sizeof(unsigned char) * h * w * 3);
  unsigned int i, j, k, l;

  for (i = 0; i < h; i++) {
    k = 3 * w * i;
    for (j = 0; j < w; j++) {
      l = 3 * j + k;
      img1D[l] = imgR[i][j];
      img1D[l + 1] = imgG[i][j];
      img1D[l + 2] = imgB[i][j];
    }
  }
  return (img1D);
}

// Eliberarea memoriei matricei
void freeImage2D(unsigned char ** mat) {
  free(mat[0]);
  free(mat);
}

int main(int argc, char * argv[]) {

  if (argc < 3) {
    printf("Wrong number of arguments. The use of the program is ./resize <image_name> <factor>\n");
    return 0;
  }

  string inPath("../../input/");
  string inFile(argv[1]);
  inPath.append(inFile);
  char outFile[50];
  float factor = atof(argv[2]);
  int cores = atoi(argv[3]);

  unsigned char * image1D;
  unsigned int w, h, bitdepth, bitsXpixel, channels, isGrey, haveAlpha;

  // Citirea imaginii
  image1D = readPNG(inPath.data(), w, h, bitdepth, bitsXpixel, channels, isGrey, haveAlpha);

  // Obtinerea fiecarui canal de culoare
  unsigned char ** imgR = getChannel2D(image1D, h, w, 0);
  unsigned char ** imgG = getChannel2D(image1D, h, w, 1);
  unsigned char ** imgB = getChannel2D(image1D, h, w, 2);

  unsigned char ** img2gR = createMatrix((int)(factor * h), (int)(factor * w));
  unsigned char ** img2gG = createMatrix((int)(factor * h), (int)(factor * w));
  unsigned char ** img2gB = createMatrix((int)(factor * h), (int)(factor * w));

  int id, r;
  void * status;
  pthread_t threads[cores];
  pthread_struct arguments[cores];
  pthread_mutex_t mutex_R, mutex_G, mutex_B;
  if (pthread_mutex_init( & mutex_R, NULL) != 0) {
    printf("Error to initialize mutex");
    return 1;
  }

  if (pthread_mutex_init( & mutex_G, NULL) != 0) {
    printf("Error to initialize mutex");
    return 1;
  }
  if (pthread_mutex_init( & mutex_B, NULL) != 0) {
    printf("Error to initialize mutex");
    return 1;
  }

  for (id = 0; id < cores; id++) {
    arguments[id].id = id;
    arguments[id].cores = cores;
    arguments[id].factor = factor;
    arguments[id].h = h;
    arguments[id].w = w;
    arguments[id].imgR = imgR;
    arguments[id].imgG = imgG;
    arguments[id].imgB = imgB;

    arguments[id].img2gR = img2gR;
    arguments[id].img2gG = img2gG;
    arguments[id].img2gB = img2gB;

    arguments[id].mutex_R = mutex_R;
    arguments[id].mutex_B = mutex_B;
    arguments[id].mutex_G = mutex_G;

    r = pthread_create( & threads[id], NULL, f, & arguments[id]);

    if (r) {
      printf("Eroare la crearea thread-ului %ld\n", id);
      exit(-1);
    }
  }

  for (id = 0; id < cores; id++) {
    r = pthread_join(threads[id], & status);

    if (r) {
      printf("Eroare la asteptarea thread-ului %ld\n", id);
      exit(-1);
    }
  }

  // Fiecare canal 2D este convertit într-o matrice 1D separată, astfel încât să poată fi salvate ca PNG
  unsigned char * imgR1 = convert2Dto1D(img2gR, factor * h, factor * w);
  unsigned char * imgG1 = convert2Dto1D(img2gG, factor * h, factor * w);
  unsigned char * imgB1 = convert2Dto1D(img2gB, factor * h, factor * w);

  pthread_mutex_destroy( & mutex_R);
  pthread_mutex_destroy( & mutex_G);
  pthread_mutex_destroy( & mutex_B);

  // Combina toate cele trei canale intr-o singura matrice 1D pentru a forma o imagine color
  unsigned char * imgRGB1 = convert2Dto1D(img2gR, img2gG, img2gB, (int)(factor * h), int(factor * w));

  // Salvarea imaginii
  sprintf(outFile, "%s%.2fx%s", "../../output/pthread/", factor, inFile.data());
  savePNG(outFile, imgRGB1, floor(factor * w), floor(factor * h));

  // Eliberare memorie
  free(image1D);
  freeImage2D(imgR);
  freeImage2D(imgG);
  freeImage2D(imgB);
  freeImage2D(img2gR);
  freeImage2D(img2gG);
  freeImage2D(img2gB);
  free(imgR1);
  free(imgG1);
  free(imgB1);
  free(imgRGB1);

  return (0);
}