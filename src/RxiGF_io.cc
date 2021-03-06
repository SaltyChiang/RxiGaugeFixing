#include "../include/RxiGF_lattice.h"
#include "../include/RxiGF_io.h"

namespace RxiGF
{

bool swapEdian = false;
size_t precision = sizeof(PREC);

PREC SwapOrder(PREC *buffer)
{
  void *bufferOut = malloc(precision);
  if (precision == 4)
  {
    long bufferInt = *(long *)buffer;
    *(long *)bufferOut = ((bufferInt & 0x000000FF) << 24) |
                         ((bufferInt & 0x0000FF00) << 8) |
                         ((bufferInt & 0x00FF0000) >> 8) |
                         ((bufferInt & 0xFF000000) >> 24);
  }
  else if (precision == 8)
  {
    long long bufferInt = *(long long *)buffer;
    *(long long *)bufferOut = ((bufferInt & 0x00000000000000FF) << 56) |
                              ((bufferInt & 0x000000000000FF00) << 40) |
                              ((bufferInt & 0x0000000000FF0000) << 24) |
                              ((bufferInt & 0x00000000FF000000) << 8) |
                              ((bufferInt & 0x000000FF00000000) >> 8) |
                              ((bufferInt & 0x0000FF0000000000) >> 24) |
                              ((bufferInt & 0x00FF000000000000) >> 40) |
                              ((bufferInt & 0xFF00000000000000) >> 56);
  }
  else
  {
    printf("Data precision not supported. \nPlease set single or double precision floating number to calculate.");
    exit(-1);
  }

  return *(PREC *)bufferOut;
}

int ReadConf(ColorMatrix *gf, char *fn)
{
#ifdef _WIN32
  FILE *fp;
  fopen_s(&fp, fn, "rb");
#else
  FILE *fp = fopen(fn, "rb");
#endif
  PREC *ptr = (PREC *)gf;
  fread((void *)ptr, 2 * precision, VOL * Nd * Nc * Nc, fp);

  if (swapEdian)
#pragma omp parallel for
    for (int i = 0; i < VOL * Nd * Nc * Nc * 2; i++)
      ptr[i] = SwapOrder(&(ptr[i]));

  fclose(fp);
  return 1;
}

int WriteConf(ColorMatrix *gf, char *fn)
{
#ifdef _WIN32
  FILE *fp;
  fopen_s(&fp, fn, "wb+");
#else
  FILE *fp = fopen(fn, "wb+");
#endif
  void *ptr = (void *)gf;
  fwrite(ptr, sizeof(Complex), VOL * Nd * Nc * Nc, fp);

  fclose(fp);
  return 1;
}

int WriteLambda(ColorMatrix *lf, char *fn)
{
#ifdef _WIN32
  FILE *fp;
  fopen_s(&fp, fn, "wb+");
#else
  FILE *fp = fopen(fn, "wb+");
#endif
  void *ptr = (void *)lf;
  fwrite(ptr, sizeof(Complex), VOL * Nd * Nc * Nc, fp);

  fclose(fp);
  return 1;
}

int WriteRandom(double *rd, char *fn)
{
#ifdef _WIN32
  FILE *fp;
  fopen_s(&fp, fn, "wb+");
#else
  FILE *fp = fopen(fn, "wb+");
#endif
  void *ptr = (void *)rd;
  fwrite(ptr, sizeof(double), VOL * (Nc * Nc - 1), fp);

  fclose(fp);
  return 1;
}

} // namespace RxiGF