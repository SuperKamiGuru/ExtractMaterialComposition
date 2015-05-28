#ifndef ISOTOPEMASS_HH
#define ISOTOPEMASS_HH

#include <string>
#include <iostream>

using namespace std;

class IsotopeMass
{
    public:
        IsotopeMass();
        virtual ~IsotopeMass();
        static void ClearStore();
        static void SetIsotopeMass();
        static double GetIsotopeMass(int Z,int A)
        {
            int iso = A-elemBaseA[Z];
            if((iso<=elemNumIso[Z])&&(iso>=0))
            {
                return isotopeMass[Z][iso];
            }
            else
            {
                cout << "\nError: isotope " << Z << " " << A << " is beyond the scope of this container" << endl;
                return 0.;
            }

        }
        static double **isotopeMass;
        static int *elemNumIso;
        static int *elemBaseA;
    protected:
    private:

};
#endif // ISOTOPEMASS_HH
