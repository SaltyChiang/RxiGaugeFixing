#pragma once

template <typename T>
void printMatrix(T &matrix)
{
  for (int i = 0; i < Nc; i++)
  {
    for (int j = 0; j < Nc; j++)
      printf("(%10.6f,%10.6f) ", matrix.elem(i, j).real(), matrix.elem(i, j).imag());
    printf("\n");
  }
}