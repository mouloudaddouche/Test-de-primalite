#include "main.h"

//permet d'affichier compose si p !=0 sinon elle affiche premier
void verification_primalite(int p){

	if (p) printf("compose\n\n");
	else printf("premier \n\n");
	
}

//Fonction Qui Permet de calculer pow(n,p) ⁼ n^p : le resultat du calcul est stocké dans la variable résultat
// elle utilisé seulement par la decomposition ( trouver le s et t dans miller rabbin)
void Puissance(mpz_t resultat,mpz_t n, mpz_t p){
	// Declaration : indice de la boucle 
	mpz_t indice;
	// Un = variable qui contient l'entier 1
	mpz_t Un; 
	// Initialisation + Affectation
	mpz_init_set_str(indice, "0", 10);
	mpz_init_set_str(Un, "1", 10); 
	mpz_init_set_str(resultat, "1", 10);
	while(mpz_cmp(p,indice)>0){
	// Multiplication : Resultat = Resultat * n
	mpz_mul(resultat,resultat,n); 
	// Incrementation
	mpz_add(indice,indice,Un); 
}
	// Liberation de la memoire allouer
	mpz_clear(indice); 
	mpz_clear(Un); 
}

// Fonction qui permet de decomposer n-1 en 2^s * t
void Decomposition (mpz_t t,mpz_t s,mpz_t n){
	mpz_t q,r,Zero,Un,Deux,indice,puissance;
	mpz_init_set_str(Zero, "0", 10);
	mpz_init_set_str(Un, "1", 10);
	mpz_init_set_str(Deux, "2", 10);
	mpz_init_set_str(indice, "1", 10);
	mpz_init(r);
	mpz_init(q);
	mpz_init(puissance);
	mpz_sub(n,n,Un);
	while(mpz_cmp(r,Zero)==0){
		Puissance(puissance,Deux,indice);
		mpz_cdiv_qr(q,r,n,puissance);
		if(mpz_cmp(r,Zero)>=0){
			mpz_set(s,indice);
			mpz_set(t,q);
		}
			mpz_add(indice,indice,Un);
	}
	mpz_clear(Zero);
	mpz_clear(Un);
	mpz_clear(Deux);
	mpz_clear(indice);
	mpz_clear(puissance);
	mpz_clear(q);
	mpz_clear(r);
}

int TestFermat(mpz_t n,int k){
	gmp_randstate_t  gmpRandState; // Random Generator
	gmp_randinit_mt(gmpRandState);
	gmp_randseed_ui(gmpRandState, time(NULL));
	mpz_t rand,tmp,Un,res;
	mpz_init(rand); 
	mpz_init(tmp);
	mpz_init(res);
	
	mpz_init_set_str(Un, "1", 10);
	for(int i=1;i<=k;i++){
	printf(" iteration = %d \n",i);
	do{
	mpz_set(tmp,n);
	mpz_sub(tmp,tmp,Un);
	mpz_urandomm(rand, gmpRandState, tmp);
	}while(mpz_cmp(rand,Un)<=0 || mpz_cmp(rand,tmp)==0);
	squareAndMultiply(res,rand,n,tmp);
	if(mpz_cmp(res,Un)!=0) {
				mpz_clear(Un);
				mpz_clear(rand);
				mpz_clear(res);
				mpz_clear(tmp);
				return 1;
			}
		}
		mpz_clear(Un);
		mpz_clear(rand);
		mpz_clear(res);
		mpz_clear(tmp);
	return 0;
}


// retourne 1 si c'est composé, 0 sinon
int TestMillerRabin(mpz_t n,int k){
	gmp_randstate_t  gmpRandState; // Random Generator
	gmp_randinit_mt(gmpRandState);
	gmp_randseed_ui(gmpRandState, time(NULL));
	mpz_t Un,Zero,Deux,s,t,y,a,N_Moins_Un,j,temporaire;
	mpz_init(s);
	mpz_init(a);
	mpz_init(t);
	mpz_init(y);
	mpz_init(N_Moins_Un);
	mpz_init(temporaire);
	mpz_init_set_str(Un, "1", 10);
	mpz_init_set_str(Zero, "0", 10);
	mpz_init_set_str(Deux, "2", 10);
	mpz_set(temporaire,n);
	mpz_set(N_Moins_Un,n);
	mpz_sub(N_Moins_Un,N_Moins_Un,Un);
	int bool;
	if ((mpz_cmp(n,Un)==0) || (mpz_cmp(n,Deux)==0)) return 1; // cas où n=1 ou n=2
	// Decomposition de n-1
	Decomposition(t,s,temporaire); 
	for(int i=1;i<=k;i++){
		printf(" iteration : %d \n",i);
		bool=0;
		mpz_init_set_str(j, "1", 10);
		do{
			mpz_urandomm(a, gmpRandState,n);  // Generer un nombre aleatoire
		}while(mpz_cmp(a,Zero)<=0 || mpz_cmp(a,n)==0);
	squareAndMultiply(y,a,n,t); // Calcule y = a^t mod n
	mpz_sub(temporaire,s,Un); // temporaire = s-1;
	if(mpz_cmp(y,Un)!=0 && mpz_cmp(y,N_Moins_Un)!=0){
		while(mpz_cmp(j,temporaire)<=0){
			squareAndMultiply(y,y,n,Deux); // Calcule y = y^2 mod n
			if(mpz_cmp(y,Un)==0) { // Cas compose
				// Liberation de la memoire
				mpz_clear(Un);
				mpz_clear(Zero);
				mpz_clear(Deux);
				mpz_clear(s);
				mpz_clear(t);
				mpz_clear(y);
				mpz_clear(a);
				mpz_clear(N_Moins_Un);
				mpz_clear(j);
				mpz_clear(temporaire);
				return 1;
			}
			if(mpz_cmp(y,N_Moins_Un)==0) { // On arrete la boucle j, en continue avec i
				bool=1;
				break;
			} 
			mpz_add(j,j,Un); // Incrementation de j ( j++ );
		}
		if (bool== 0 ){ // Cas composé
			// Liberation de la memoire
			mpz_clear(Un);
			mpz_clear(Zero);
			mpz_clear(Deux);
			mpz_clear(s);
			mpz_clear(t);
			mpz_clear(y);
			mpz_clear(a);
			mpz_clear(N_Moins_Un);
			mpz_clear(j);
			mpz_clear(temporaire);
			return 1;
		}	
	}
  }
  // Cas Premier
  // Liberation de la memoire
  mpz_clear(Un);
  mpz_clear(Zero);
  mpz_clear(Deux);
  mpz_clear(s);
  mpz_clear(t);
  mpz_clear(y);
  mpz_clear(a);
  mpz_clear(N_Moins_Un);
  mpz_clear(j);
  mpz_clear(temporaire);
return 0;
}

// Fonction qui calcule exponentiation modulaire rapide
void squareAndMultiply( mpz_t res, mpz_t a,mpz_t n,mpz_t H){	
	mpz_t Zero,Deux,Un,r,Dix,temp,mod;
	mpz_init_set_str(Zero, "0", 10);
	mpz_init_set_str(Deux, "2", 10);
	mpz_init_set_str(Un, "1", 10);
	mpz_init_set_str(r, "1", 10);
	mpz_init_set_str(Dix, "10", 10);
	mpz_init (temp);
	mpz_init (mod);
	mpz_set(temp,H);
	while(mpz_cmp(temp,Zero)>0){
			mpz_mod(mod,temp,Deux);
			if(mpz_cmp(mod,Zero)!=0) {
			mpz_mul(r,r,a);
			mpz_mod(r,r,n);}
			mpz_mod(a,a,n);
			mpz_mul(a,a,a);
			mpz_mod(a,a,n);
		    mpz_fdiv_q(temp,temp,Deux);
}
	mpz_set(res,r);
}


// La fonction Main qui contiendra notre Menu Principale de l'application avec les differentes options
int main() {
srand(time(NULL));
mpz_t n,Zero,Un,a,b;
mpz_init(n);
mpz_init(a);
mpz_init(b);
int k=0;
int res=0;
int choix_menu_1;
int choix_menu_2;
int choix_menu_4;
int menu_1;
int menu_2;
int menu_3;
int menu_4;
mpz_init_set_str(Zero, "0", 10);
mpz_init_set_str(Un, "1", 10);
Accueil :
system ("clear");
printf(" \n ********* Bienvenue Dans Application : Test De Primalite ********* \n\n");
choix_menu_1=0;
choix_menu_2=0;
choix_menu_4=0;
menu_1=0;
menu_2=0;
menu_3=0;
menu_4=0;
printf(" veuillez choisir le test \n");
printf(" 1) Test de Fermat \n");
printf(" 2) Test de Miller-Rabin \n");
do{
 if (menu_1==0){
	scanf("%d",&choix_menu_1);
	menu_1 =1;}
	else{
		printf(" Veuillez faire rentrer un choix correcte ( 1 ou 2 ) \n");
		scanf("%d",&choix_menu_1);}
	}while((choix_menu_1 != 1)  && (choix_menu_1 != 2));
	
	printf(" Veuillez faire rentrer l'entier : n \n");
	printf(" 1) Format basique \n");
	printf(" 2) Format a^b - 1 \n");

   do{
	if (menu_2==0){
	scanf("%d",&choix_menu_2);
	menu_2 =1;}
	else{
		printf(" Veuillez faire rentrer un choix correcteaa ( 1 ou 2 ) \n");
		scanf("%d",&choix_menu_2);}
	}while((choix_menu_2 != 1)  && (choix_menu_2 != 2));

	if(choix_menu_2==1){
	do{
	if (menu_2==0){
	gmp_scanf("%Zd",&n);
	menu_2=1;
	}
	else{
	printf(" Veuillez faire rentrer un entier correcte ( n > 0) \n");
	gmp_scanf("%Zd",&n);
    }
  }while(mpz_cmp(n,Zero)<=0);
}

if(choix_menu_2==2){
	printf(" Veuillez faire rentrer a  (a > 1)\n");
	do{
		gmp_scanf("%Zd",&a);
	}while(mpz_cmp(a,Un) <=0);
	
	printf(" Veuillez faire rentrer b ( b > 1 )\n");
	do{
		gmp_scanf("%Zd",&b);
	}while(mpz_cmp(b,Un) <=0);
	Puissance(n,a,b);
	mpz_sub(n,n,Un);
	gmp_printf(" n = %Zd \n",n);
	}

printf(" Veuillez faire rentrer le nombre de repetitions : k \n");
do{
	if (menu_3==0){
	gmp_scanf("%d",&k);
	menu_3=1;
	}
else{
	printf(" Veuillez faire rentrer un entier correcte ( k > 0) \n");
	gmp_scanf("%d",&k);
}
}while(k<=0);

Menu_choix :
if(choix_menu_1 ==1){ // Test Fermat
	res=TestFermat(n,k);
}

if ( choix_menu_1 ==2){ // Test Miller-Rabbin
	res=TestMillerRabin(n,k);
	
}
	printf(" Le nombre saisie est un entier ");
	verification_primalite(res);
	
	printf(" Voulez vous ? : \n");
	printf(" 1) Retourner a l'accueil \n");
	printf(" 2) Utiliser la meme configuration pour effectuer le test ");
	if (choix_menu_1==1) {
		printf("Miller-Rabbin \n");
		choix_menu_1 =2;
		}
	else {
		printf("Fermat\n");
		choix_menu_1=1;
		}
	printf(" 3) Quiiter \n");
	
	do{
	if (menu_4==0){
	gmp_scanf("%d",&choix_menu_4);
	menu_4=1;
	}
else{
	printf(" Veuillez faire rentrer un choix correcte ( 1 , 2 ou 3) \n");
	gmp_scanf("%d",&choix_menu_4);
}
}while((choix_menu_4 != 1) && (choix_menu_4 != 2) && (choix_menu_4 != 3));

if (choix_menu_4 == 1) goto Accueil;
if (choix_menu_4 ==2 ) {
	system ("clear");
	menu_4=0;
	goto Menu_choix;}
if(choix_menu_4 ==3) {
	system ("clear");
	mpz_clear(n);
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(Zero);
	mpz_clear(Un);
	printf("             \n\n\n\n *********** Au Revoir *********** \n\n\n\n");
	return 0;}
	
}

