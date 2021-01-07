#ifndef _MATRICE
#define _MATRICE
#include <iostream>
#include <ctime>
#include <cmath>
#include "matrice.h"
#include <math.h>
#include <vector>
#include <algorithm> 

const double SMALL = 1.0E-30;


using namespace std;

template <class T>
class matrice{
		
		public: 
		
			// costruttori
			matrice();
			matrice(int,int);
			
			// distruttore
			~matrice(); 		
			
			 // costruttore di copia
			matrice(const matrice<T>&); 
			
			// sovraccarico operatore = per l'assegnamento
			matrice<T> & operator=(matrice<T> const &); 
			
			// sovraccarico operatore * per il prodotto
			matrice<T> operator*( matrice<T>&); 

			//  visitatori  	     
			int getRighe (void);	           
			int getColonne(void);
			T getDet();
			//tipoelem getDet(void);
			
			// assegnatori              
			void setRighe(int);                
			void setColonne(int);			
			
			// ritorna il valora della matrice nella posizione indicata
			//tipoelem leggimatrice(int,int);  
			T leggimatrice(int,int);  
			
			 // scrive nella matrice il valore indicato nella posizione indicata
			void scrivimatrice(int,int,T);
			
			// ritorna lo scalare della matrice
			matrice<T> scalare(int); 		     
			
			// ritorna la trasposta della matrice
			matrice<T> trasposta();		     
			
			// ritorna il prodotto tra la due matrici
			matrice<T> prodotto(matrice&);	     
			
			//  calcola il determinante della matrice con il metodo di Laplace
			T Det(int);		// deprecato, non usarlo
			
			
			//  calcola il determinante con il metodo di Gauss
			double gaussDet(void);	// da rifare
			
		        //  stampa la matrice  
			void stampa();
			
			//  scrive la matrice con valori casuali in un range
			void asRand(int,int);  
			
			
			
			             
			
		private:
			int righe;
			int colonne;
			//tipoelem det;
			//tipoelem **elementi;
			double det;
			T **elementi;
};



// costruttore
template <class T>
matrice<T>::matrice(int c,int r){

	colonne = c;

	righe = r;

	// allocazione dinamica della matrice

	int i;

	elementi = new T* [righe];

	for (i = 0; i < righe; i++)

		elementi[i] = new T[colonne];

}

// distruttore
template <class T>
matrice<T>::~matrice()
{
for (int i=0; i<righe; i++)
 
        delete[] elementi[i];
    
    delete[] elementi;
}

// costruttore di copia
template <class T>
matrice<T>::matrice(const matrice<T> &m){

	righe = m.righe;

	colonne = m.colonne;

  int i,j;

  elementi = new T* [righe]; 

  for (i=0; i < righe; i++)

    elementi[i] = new T[colonne];

  for (i=0;i<righe;i++)

		for (j=0;j<colonne;j++)

			elementi[i][j]=m.elementi[i][j];
	
}
// sovraccarico operatore =
template <class T>
matrice<T> &matrice<T>::operator=(const matrice<T> &m){

	// evita gli auto assegnamenti

	if (this == &m) return *this;

	else {

		int i,j;

		if (colonne != m.colonne || righe != m.righe){

			this->~matrice();

			colonne = m.colonne;

			righe = m.righe;

			elementi = new T* [righe]; 

			for (i=0; i < righe; i++)

				elementi[i] = new T[colonne];

		}

		for (i=0;i<righe;i++)

			for (j=0;j<colonne;j++)

				elementi[i][j] = m.elementi[i][j];

	}

	return *this;
}

template <class T>
matrice<T> matrice<T>::operator*( matrice<T> &B){


	matrice C(getRighe(),B.getColonne());

	int i, j, k;

	T s;

	for (i=0;i<getRighe();i++) 

	     for(j=0;j<B.getColonne();j++) {

	            C.scrivimatrice(i,j,0);

        	 for(k=0;k<B.getRighe();k++){  

    	           s = C.leggimatrice(i,j) + leggimatrice(i,k) * B.leggimatrice(k,j);

    	           C.scrivimatrice(i,j,s); 

    	         }

         }
        
        return C;	

}

template <class T>
void matrice<T>::setRighe(int r)
{
    	righe=r;
}


template <class T>
void matrice<T>::setColonne(int c)
{
    	colonne=c;
}

template <class T>
int matrice<T>::getRighe()
{
	return righe;
}

template <class T>
int matrice<T>::getColonne()
{
	return colonne;
}

template <class T>
T matrice<T>::getDet()
{

	return det;

}

template <class T>
T matrice<T>::leggimatrice(int i, int j)
{

	return elementi[i][j];

}

template <class T>
void matrice<T>::scrivimatrice(int i, int j, T a)
{

	elementi[i][j] = a;

} 

template <class T>
matrice<T> matrice<T>::scalare(int k)
{

	int i, j;

	matrice B(getColonne(),getRighe());

	for (i = 0;i < B.getRighe(); i++)

		for (j = 0;j < B.getColonne(); j++)

			B.scrivimatrice(i,j,leggimatrice(i,j)*k);

	return B;  
}

template <class T>
matrice<T> matrice<T>::trasposta(){

	int i, j;

	matrice B(getRighe(),getColonne());

	for (i = 0;i < righe; i++)

		for (j = 0;j < colonne; j++)

			B.scrivimatrice(j,i,leggimatrice(i,j));		

	return B;

}

template <class T>
matrice<T> matrice<T>::prodotto(matrice<T>& B){	

	matrice C(getRighe(),B.getColonne());

	int i, j, k;

	T s;

	for (i=0;i<getRighe();i++) 

	     for(j=0;j<B.getColonne();j++) {

	            C.scrivimatrice(i,j,0);

        	 for(k=0;k<B.getRighe();k++){  

    	           s = C.leggimatrice(i,j) + leggimatrice(i,k) * B.leggimatrice(k,j);

    	           C.scrivimatrice(i,j,s); 

    	         }

    	     }	


	return C;

}

template <class T>
void matrice<T>::stampa(){

for (int i=0;i<righe;i++){

		cout << "\n";

		for (int j=0;j<colonne;j++){	

			cout << leggimatrice(i,j) << "  " ;

		}

	}	

	cout << "\n";


}


template <class T>
T matrice<T>::Det(int n) {	// deprecato
   det = 0;
   matrice submatrix(n,n);
   if (n == 2)
   return ((leggimatrice(0,0)*leggimatrice(1,1))-(leggimatrice(1,0)*leggimatrice(0,1)));
   else {
      for (int x = 0; x < n; x++) {
         int subi = 0;
         for (int i = 1; i < n; i++) {
            int subj = 0;
            for (int j = 0; j < n; j++) {
               if (j == x)
               continue;
               submatrix.scrivimatrice(subi,subj,leggimatrice(i,j));
               subj++;
            }
            subi++;
         }
         det = det + (pow(-1, x) * leggimatrice(0,x) * submatrix.Det(n - 1 ));
      }
   }
   return det;
} 


template <class T>
void matrice<T>::asRand(int start,int range){
      	
      int i, j;

      if (typeid(T) == typeid(double)){		//TODO : non funziona il limite massimo e minimo
          	
      	double f;

	      for (i=0;i<getRighe();i++){		
		
			for (j=0;j<getColonne();j++){
			
				f = (double)rand() / RAND_MAX;
				
				scrivimatrice(i,j,range + f * (range - start));
			
			}
			
		}
    
      }


       else{
    

	      for (i=0;i<getRighe();i++){		
		
			for (j=0;j<getColonne();j++){
				
				scrivimatrice(i,j,start + rand() % range);
			
			}
			
		}

        }

}

template <class T>
double matrice<T>::gaussDet(void){	// da rifare

   int i, j, k, r;
   double maxA, val, pivot, multiple;
   det = 1;

   for ( i = 0; i < righe - 1; i++ )
   {
      r = i;
      maxA = abs( leggimatrice(i,i) );
      for ( k = i + 1; k < righe; k++ )
      {
         val = abs( leggimatrice(k,i) );
         if ( val > maxA )
         {
            r = k;
            maxA = val;
         }
      }
      if ( r != i )
      {
         for ( j = i; j < righe; j++ ) swap( elementi[i][j], elementi[r][j] );
         det = -det;
      }
      
      

      pivot = leggimatrice(i,i);
      if ( abs( pivot ) < SMALL ) return 0.0;              // matrice singolare

      for ( r = i + 1; r < righe; r++ )                    
      {
         multiple = leggimatrice(r,i) / pivot;                
         for ( j = i; j < righe; j++ ) scrivimatrice(r,j,(leggimatrice(r,j) - (multiple * leggimatrice(i,j))));  
      }
      det *= pivot;                                        
   }

   det *= leggimatrice(righe-1,righe-1); 

   return det;

}



#endif

