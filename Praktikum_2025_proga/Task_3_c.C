#include <iostream>
#include <fstream>
#include <math.h>

using namespace ::std;

struct Point{
  double x;
  double y;
};
struct Triangle{
  struct Point A;
  struct Point B;
  struct Point C;

};
void fillTriangle(double Ax, double Ay, double Bx, double By, double Cx, double Cy, struct Triangle *triangle_to_fill)
{
  triangle_to_fill->A.x = Ax;
  triangle_to_fill->A.y = Ay;
  triangle_to_fill->B.x = Bx;
  triangle_to_fill->B.y = By;
  triangle_to_fill->C.x = Cx;
  triangle_to_fill->C.y = Cy;
}
void printTriangle(struct Triangle t)
{
  printf("A(%.2f, %.2f)\n", t.A.x, t.A.y);
  printf("B(%.2f, %.2f)\n", t.B.x, t.B.y);
  printf("C(%.2f, %.2f)\n", t.C.x, t.C.y);
}
double calculateArea(struct Triangle t) {
    double area = 0.5 * fabs(
        (t.A.x * (t.B.y - t.C.y) +
         t.B.x * (t.C.y - t.A.y) +
         t.C.x * (t.A.y - t.B.y))
    );
    return area;
}
struct Vectr{
  double x;
  double y;
};
void fillVectr(struct Vectr *vectr, struct Point A, struct Point B) //vector AB
{
  vectr->x = B.x - A.x;
  vectr->y = B.y - A.y;
}
void rotateVectorToPi(struct Vectr vectr_to_rotate, struct Vectr *rotated_vectr)
{
  rotated_vectr->x = vectr_to_rotate.y;
  rotated_vectr->y = -vectr_to_rotate.x;
}
double scalarProduct(struct Vectr vector_1,struct Vectr vector_2)
{
  return vector_1.x*vector_2.x + vector_1.y*vector_2.y;
}

bool isPoinrInsideTriangle(struct Triangle triangle, struct Point point)
{
  struct Vectr AB;
  struct Vectr BC;
  struct Vectr CA;

  fillVectr(&AB, triangle.A,triangle.B);
  fillVectr(&BC, triangle.B,triangle.C);
  fillVectr(&CA, triangle.C,triangle.A);

  //now check is the point by one side with third point of triangle
  //checking point and third point need to have same sign of projection to transporant line
  
  struct Vectr AB_T;
  struct Vectr BC_T;
  struct Vectr CA_T;

  rotateVectorToPi(AB,&AB_T);
  rotateVectorToPi(BC,&BC_T);
  rotateVectorToPi(CA,&CA_T);

  struct Vectr A_point;
  struct Vectr B_point;
  struct Vectr C_point;

  fillVectr(&A_point,triangle.A,point);
  fillVectr(&B_point,triangle.B,point);
  fillVectr(&C_point,triangle.C,point);

  bool is_Apoint_AC_AB_T = scalarProduct(A_point,AB_T)*scalarProduct(CA,AB_T)<0.;
  bool is_Bpoint_BA_BC_T = scalarProduct(B_point,BC_T)*scalarProduct(AB,BC_T)<0.;
  bool is_Cpoint_CB_CA_T = scalarProduct(C_point,CA_T)*scalarProduct(BC,CA_T)<0.;

  return is_Apoint_AC_AB_T*is_Bpoint_BA_BC_T*is_Cpoint_CB_CA_T;
}
bool isZeroInsideTriangle(struct Triangle triangle)
{
  struct Point Zero_Point;
  Zero_Point.x=0.;
  Zero_Point.y=0.;
  bool a = isPoinrInsideTriangle(triangle,Zero_Point);
  return a;
}

int compareArea(const void* a, const void* b) {
    struct Triangle *triA = (struct Triangle *)a;
    struct Triangle *triB = (struct Triangle *)b;
    
    double areaA = calculateArea(*triA);
    double areaB = calculateArea(*triB);
    
    if (areaA < areaB) return -1;
    if (areaA > areaB) return 1;
    return 0;
}

void Task_3_c()
{
  FILE *file = fopen("Praktikum_2025_proga/triangles.txt","r");
  if (!file)
  {
    printf("Ошибка открытия файла\n");
    return ;
  }

  const int N_arr_size = 1000;
  struct Triangle Arr_Triangles[N_arr_size];

  int count = 0;
  const int stop_check = 10000;

  while (fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf",
                &Arr_Triangles[count].A.x, &Arr_Triangles[count].A.y,
                &Arr_Triangles[count].B.x, &Arr_Triangles[count].B.y,
                &Arr_Triangles[count].C.x, &Arr_Triangles[count].C.y) == 6 && count<stop_check)
  {
    count++;
  }

  qsort(Arr_Triangles, count, sizeof(Triangle), compareArea);

  double Area_min,Area_max;

  Area_min = calculateArea(Arr_Triangles[0]);
  Area_max = calculateArea(Arr_Triangles[N_arr_size-1]);

  printf("Minimum area %f of triangle:\n", calculateArea(Arr_Triangles[0]));
  printTriangle(Arr_Triangles[0]);

  printf("Maximum area %f of triangle:\n", calculateArea(Arr_Triangles[N_arr_size-1]));
  printTriangle(Arr_Triangles[N_arr_size-1]);

  int N_of_zero_inside=0;
  for(int i=0;i<N_arr_size;i++)
  {
    if(isZeroInsideTriangle(Arr_Triangles[i]))
    {
      N_of_zero_inside++;
    }
  }

  printf("\nNumber of triangles with zero inside: %d \n\n", N_of_zero_inside);

  fclose(file);
}