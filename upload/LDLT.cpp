/*
 *  RATEFIARISON 
 *  Harivony Lalatiana
 * 
 *  Factorisation LDLt
 * 
 */

#include <iostream>     // bibliothèque standard d'entrE et de sortie
#include <fstream>      // bibliothèque d'entrE et de sortie pour les fichier
#include <string>       // bibliotheque de manipulation de chaine de charactere
#include <vector>
#include <iomanip>

using namespace std;    // utiliser l'espace de nom standard std


int dim(10);            // definir la taille du systeme 

/*
    Allocation/Affichage de matrice et de vecteur
*/

float **newMat(int rows, int cols);     // creer une matrice rows x cols
float *newVec(int n);                   // creer un vecteur de taille n
void displayMat(float** mat);           // afficher une matrice de mat
void displayVec (float* v);             // afficher un vecteur 

/*
 *
 *  CLASS: LDLt
 * 
 * 
*/

class LDLT{
// Methods
public:
    LDLT();
    ~LDLT();
    void loadAb();
    void computeP_i();
    void computeAp_nDiag();
    void computeAP();
    void computeX();

    // getter
    float** getA(){return A;}
    float* getb(){return b;}
    float** getD(){return D;}
    float* getX(){return x;}
    float getAp(int i, int j);
    void setAp(int i, int j, float val);

    void solveTriangInf();
    void solveDiag();
	void solveTriangSup();

    void profil(float** mat, string name);
    void getProfil(float* mat, string name);
    void solve();
    void gaussSolve();
    vector<float> AP = {};
	vector<int> l = {};
	vector<int> p = {};
    vector<int> nDiag;

// Attributs
private:
    float** A;
    float** L;
    float** D;
    float* x;    // la solution

    float* b;    // second membre
    float* p_i;  // indice des premieres element non null
};

/*
 *
 *  MAIN
 * 
 * 
*/

int main (){
    // initialiser une instance de LDLt 
    LDLT solver;

    // Interface Utilisateur
    cout << "____________________________________" << endl;
    cout << "\n  Factorisation LDLt et Resolution " << endl;
    cout << "____________________________________" << endl << endl;

    cout << "Factoriser une matrice symetric positive A\n tq: A = L.D.Lt"<< endl<< endl;
    // factorisation et resolution
    solver.solve();

    // succes de l'execution
    return 0;
}

/***** MATRICE - VECTEUR ******/

float **newMat(int rows, int cols){ //allocation dynamique de tableau rows*cols
    float **mat(NULL);
    mat = new float*[rows];

    for (int i=0;i<rows;i++)        
        mat[i] = new float[cols];
    return mat;
}

float *newVec(int n){               //alllouer un vecteur de dim n
    float* b;
    b  = new float [n];
    return b;
}

void displayMat(float** mat){    //afficher la matrice
    for (int i=0;i<dim;i++){
        for(int j=0; j<dim; j++)
            cout << mat[i][j] << " ";
        cout << endl;
    }
}

void displayVec (float* v){         //afficher un vecteur
    for(int i=0; i<dim; i++){
        cout << "[" << v[i] << "]" << endl;
    }
}


/*
 *
 *  DEFINITION : METHODE LDLT
 * 
 * 
*/

LDLT::LDLT(){
    // Allcation des matrice et des vecteurs
    A = newMat(dim,dim);
    L = newMat(dim,dim);
    D = newMat(dim,dim);
    b = newVec(dim);
    x = newVec(dim);
}

// Chargger la matrice depuis un fichier
void LDLT::loadAb(){
    string fileName;
    ifstream inFile;
    
    //cout << "Entrez le fichier a charger : "; cin >> fileName;

    // ouverture du fichier, 
    inFile.open("matrix.txt");

    if(inFile){
        //recuperer la taille de la matrice
        inFile >> dim;

        // recuperer la partie inférieure avec le duagonal de la matrice
        for(int i=0; i<dim;i++){
            for(int j=0; j<dim;j++){
                inFile >> A[i][j];
            }
        }  
        
        // recuperer le second membre
        for(int i=0; i<dim;i++){
                inFile >> b[i];
        }
    }
    
    // Alerter l'erreur d'ouverture du fichier
    else
        cout << "erreur de lecture" << endl;
}

// charger le profil p_i de A
void LDLT::computeP_i(){
    p_i = newVec(dim);

    for(int i=0;i<dim;i++){
        bool logicNull = true;
        for(int j=0;j<=i; j++){
            if(A[i][j]==0){
                if(logicNull)
                    cout <<"  ";
                else
                    cout <<". ";
            }
            else {
                cout << "x ";
                if(logicNull)
                    p_i[i] = j;        // enregistrer la valeur de j comme etant le premmier element non null de la ligne
                logicNull = false;     // les zeros qui suive ne sont plus null logiquement
            }
        }
        cout << endl;
    }
    cout << endl;
}


// afficher le profi de la matrice mat
void LDLT::profil(float** mat, string name){
    cout << "Profil de la matrice " << name << " :"<< endl;
    for(int i=0;i<dim;i++){
        // eviter de faire un boucle de plus
        bool logicNull = true;
        for(int j=0;j<=i; j++){
            if(mat[i][j]==0){
                if(logicNull)
                    continue;
                else
                    cout <<"("<< i+1 << ";" << j+1 << "), ";
            }
            else {
                cout <<"("<< i+1 << ";" << j+1 << "), ";
                logicNull = false; // les zeros qui suive ne sont plus null logiquement
            }
        }
        cout << endl;
    }
    cout << endl;
}

void LDLT::computeAp_nDiag(){
	ifstream file("matrix.txt", ios::in);
    if (file){
        file >> dim;

        // Initialisations depuis le fichier
        int ji(0);
        int count(-1); // Une variable qui va numeroter les éléments du profil
        float elt(0);
        char car('c');

        for(int_fast16_t i=0; i<dim; i++){
            ji = 0; 
            elt = 0;

            // Prendre seulement les éléments du profil
            for(int_fast16_t j=0; j<=i; j++){
                file >> elt;
                if (elt != 0 || ji != 0){
                    ji++;
                    if(ji == 1){
                        p.push_back(j);
                        l.push_back(i-j);
                    }
                }
                if(ji){
                    count++;
                    AP.push_back(elt);
                    if( i == j){
                        nDiag.push_back(count);
                    }
                }
            }            
            // Boucler jusqu'au dernier élément de la ligne
            while(car != '\n'){
                file.get(car);
            }
            car = 'c';
        }

        for(int_fast16_t i=0; i<dim; i++){
            file >> b[i];
        }
    }
}

// determiner la valeur de L et D avec A = L.D.Lt
void LDLT::computeAP(){     
    float sum = 0;
    for (int i=0; i<dim; i++){
        for(int j=0; j<i; j++){
            sum = 0;
            for(int k=0; k<j; k++){
                if (k>=p[i] && k>= p[j])
                    sum += AP[nDiag[i]-i+k]*AP[nDiag[k]]*AP[nDiag[j]-j+k];
            }
            if(j>=p[i] && j>=p[j])
                AP[nDiag[i]-i+j] = (1/AP[nDiag[j]])*(AP[nDiag[i]-i+j] - sum);
        }
        sum = 0; 
        for (int k=0; k<i; k++){
            if (k>=p[i])
                sum += AP[nDiag[k]]*(AP[nDiag[i]-i+k]*AP[nDiag[i]-i+k]);
        }
        AP[nDiag[i]] = AP[nDiag[i]] - sum;
    }
}

void LDLT::computeX(){
    // D'abord Lx = b
    float sum(0);
    for(int i=0; i<dim; i++){
        sum = 0;
        for(int j=0; j<i; j++){
            if(j>=p[i])
                sum += AP[nDiag[i]-i+j]*x[j];
        }
        x[i] = b[i] - sum;
    }
    // Puis Dx = x
    for(int i=0; i<dim; i++){
        x[i] = (1/AP[nDiag[i]])*x[i];
    }
    // Enfin Ltx = x
    for(int i= int(dim-1); i>=0; i--){
        sum = 0;
        for(int j=i+1; j<dim; j++){
            if(i>=int(p[j]))
                sum += AP[nDiag[j]-j+i]*x[j];
        }
        x[i] = x[i]-sum;
    }
}

float LDLT::getAp(int i, int j){
	
	if(j >= p[i] && j<=i){
		return AP[nDiag[i] - i + j];
	}	
	else 
	return 0;
}

void LDLT::setAp(int i, int j, float val){
	if(j >= p[i]) AP[nDiag[i] - i + j] = val;
}

// Appeler par etape les methodes de resolutions
void LDLT::solve(){
    // 1 - charger les donner dempuis les fichier
    loadAb();
    
    // Interface
    cout << "La matrice A :" << endl;
    displayMat(A);
    
    cout << "\nLe second membre :" << endl;
    displayVec(b);

    cout << "\nProfil de A :" << endl;
    computeP_i(); // tracer le profil

    // 2 - calcule de nDiag, l et p
    computeAp_nDiag();
    // 3 - Ajouter les Valeur de AP 
    computeAP(); 
    
    // 4 - Calculen de la solution
    computeX();  

    cout << "La solution du systeme est :"<< endl;
    displayVec(x);                   // afficher le resultat
}

// liberation de memoire apres execution
LDLT::~LDLT(){
	delete[] x;
	delete[] b;
    l.clear();
    p.clear();

    AP.clear();
    nDiag.clear();
	for(int i=0; i< dim; i++){
        delete[] D[i];
        delete[] A[i];
        delete[] L[i];
    }
    delete[] D;
    delete[] A;
    delete[] L;
}