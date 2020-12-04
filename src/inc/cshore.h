#ifndef CSHORE_H
#define CSHORE_H

extern "C"
{
   void cshore(int*, long const*, int const*, int const*);

   void CShoreWrapper(long const*, int const*, int const*, int const*, int const*, int const*, int const*, int const*, int const*, int const*, int const*, int const*, int const*, int const*, double const*, double const*, double[], double[], double[], double[], double[], double[], int const*, double[], double[], double[], int*, int*, double[], double[], double[], double[]);
}
#endif // CSHORE_H
