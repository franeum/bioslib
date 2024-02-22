#include "ext.h"
#include "ext_obex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <time.h>
#include <stdint.h>
#include "mt64.h"

// strutture

typedef struct _cost
{
    double singlefit;
    t_atom value[16];
} t_cost;


typedef struct _BGA
{
    t_object ob;
    long var_size; // dimensione in bit delle variabili
    long pop_size; // dimensione della popolazione
    long n_var; // numero di variabili per individuo
    short n_integers; // numero di interi a 64 bit per rappresentare un individuo
    double lambda_in; // valore di fitness in entrata dal lambda inlet
    short matingtype; // tipo di algoritmo per l’accoppiamento
    long nbit; // numero di bit da mutare per ogni individuo
    double prob; // probabilità di mutazione
    long poplen; // lunghezza della popolazione (popsize * n_integers)
    long crosstype;

    long generationcounter;
    long refreshnum; // numero di generazioni dopo le quali fare il refresh
    long totrefresh; // numero di generazioni dopo le quali fare il refresh totale

    t_atom *population; // array per la popolazione
    t_atom *phenotype; // array per le variabili decodificate
    t_cost *individual; // array con le coppie individui/fitness

    short chromlen;

    uint64_t *mask1;
    uint64_t *mask2;

    // outlets

    void *b_lambdaout;
    void *b_out0;
    void *b_out1;
    void *b_out2;
    void *b_out3;
} t_BGA;


// prototipi di funzione

void *BGA_new(t_symbol *s, long argc, t_atom *argv);
void BGA_free(t_BGA *x);
void BGA_assist(t_BGA *x, void *b, long io, long index, char *s);
void BGA_bang(t_BGA *x);
void BGA_ft1(t_BGA *x, double fit);
void BGA_init(t_BGA *x);
void BGA_refresh(t_BGA *x, long r);
void BGA_reinit(t_BGA *x);
void BGA_totrefresh(t_BGA *x, long r);
void BGA_reinitnew(t_BGA *x);
void BGA_gen2phen(t_BGA *x);
void BGA_bubblesort(t_BGA *x);
void BGA_couples(t_BGA *x);
void BGA_crossoverSP(t_BGA *x, short mother, short father, short counter);
void BGA_crossoverDP(t_BGA *x, short mother, short father, short counter);
void BGA_crossoverU(t_BGA *x, short mother, short father, short counter);
void BGA_crossoverAN(t_BGA *x, short mother, short father, short counter);
void BGA_mutation(t_BGA *x);
void BGA_nbit(t_BGA *x, long t);
void BGA_prob(t_BGA *x, double t);
void BGA_applyMut(t_BGA *x, short index);
void BGA_repack(t_BGA *x);
void BGA_ctype(t_BGA *x, long t);

// puntatore globale alla classe

t_class *BGA_class;

//int C74_EXPORT main(void)
//void ext_main(void *r)
void ext_main(void *r)
{
    t_class *c;

    c = class_new("bios.BGA", (method)BGA_new, (method)BGA_free, (long)sizeof(t_BGA),
    0L, A_GIMME, 0);

    class_addmethod(c, (method)BGA_assist, "assist", A_CANT, 0);
    class_addmethod(c, (method)BGA_bang, "bang", 0);
    class_addmethod(c, (method)BGA_ft1, "ft1", A_FLOAT, 0);
    class_addmethod(c, (method)BGA_init, "init", 0);
    class_addmethod(c, (method)BGA_nbit, "nbit", A_LONG, 0);
    class_addmethod(c, (method)BGA_prob, "prob", A_FLOAT, 0);
    class_addmethod(c, (method)BGA_ctype, "ctype", A_LONG, 0);
    class_addmethod(c, (method)BGA_refresh, "refresh", A_LONG, 0);
    class_addmethod(c, (method)BGA_totrefresh,"totrefresh", A_LONG, 0);

    CLASS_ATTR_LONG(c, "popsize", 0, t_BGA, pop_size);
    CLASS_ATTR_LONG(c, "numvar", 0, t_BGA, n_var);
    CLASS_ATTR_LONG(c, "varsize", 0, t_BGA, var_size);

    class_register(CLASS_BOX, c); /* CLASS_NOBOX */
    BGA_class = c;

    post("∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗");
    post("−−− bios.BGA");
    post("GA Context − classical binary GA");
    post("2014 July − Francesco Bianchi");
    post("v0.3.2 alpha");
    post("∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗");

    //return 0;
}

/* inlets: 1. bang, init, e metodi vari
2. pop size
3. num vars
4. var size
5. lambda inlet */

void BGA_ft1(t_BGA *x, double fit)
{
    x->lambda_in = fit;
}


// ================================
// numero di bit per la mutazione
// ================================

void BGA_nbit(t_BGA *x, long t)
{
    x->nbit = t;
}


// =======================================
// probabilità di mutazione via messaggio
// =======================================

void BGA_prob(t_BGA *x, double t)
{
    x->prob = t;
}


// ===================
// tipo di crossover
// ===================

void BGA_ctype(t_BGA *x, long t)
{
    x->crosstype = t;
}


// ================================================
// numero di generazioni dopo cui fare il refresh
// ================================================

void BGA_refresh(t_BGA *x, long r)
{
    x->refreshnum = r;
}


// ======================================================
// numero di generazioni dopo cui fare il refresh totale
// ======================================================

void BGA_totrefresh(t_BGA *x, long r)
{
    x->totrefresh = r;
}


// ==============================================
// init: genero la popolazione iniziale (random)
// ==============================================

void BGA_init(t_BGA *x)
{
    uint64_t seed1, seed2, seed3, seed4; // semi per la generazione dei numeri random
    short i; // contatore per i cicli
    short chromlength; // dimensione di un individuo (in bit)
    short nintegers; // numero di interi per rappresentare un individuo
    short poplen; // lunghezza della popolazione
    short nmask; // quantità di 1 nella mask (da sinistra)
    short popsize;
    uint64_t mask; // bitmask per la rappresentazione dell’ultimo intero (quello incompleto)

    seed1 = rand()%0XFFFFFFFFULL;
    seed2 = rand()%0XFFFFFFFFULL;
    seed3 = rand()%0XFFFFFFFFULL;
    seed4 = rand()%0XFFFFFFFFULL;

    x->generationcounter = 0;
    popsize =x->pop_size;
    chromlength = x->n_var * x->var_size;
    nintegers = chromlength / 64;
    nmask = chromlength % 64;
    if (nmask) nintegers++;

    x->chromlen = chromlength;

    poplen = x->pop_size * nintegers;

    x->n_integers = nintegers;
    x->poplen = poplen;

    // inizializzo i semi per la generazione dei numeri casuali

    uint64_t init[4]={seed1, seed2, seed3, seed4}, length=4;
    init_by_array64(init, length);

    for (i=0; i<poplen; i++) {
        atom_setlong(&x->population[i], genrand64_int64());
    }

    if (nmask) {
        mask = UINT64_MAX << (uint64_t)(64 - nmask);
        
        for (i=1; i<=popsize; i++) {
            atom_setlong(&x->population[i * nintegers - 1], atom_getlong(&x->population[i
            * nintegers - 1]) & mask);
        }
    }

    BGA_gen2phen(x);
}


// =========================================================
// reinit: rigenerazione della popolazione dopo il refresh
// =========================================================

void BGA_reinit(t_BGA *x)
{
    uint64_t seed1, seed2, seed3, seed4; // semi per la generazione dei numeri random
    short i; // contatore per i cicli
    short chromlength; // dimensione di un individuo (in bit)
    short nintegers; // numero di interi per rappresentare un individuo
    short poplen; // lunghezza della popolazione
    short nmask; // quantità di 1 nella mask (da sinistra)
    short popsize;
    uint64_t mask; // bitmask per la rappresentazione dell’ultimo intero (quello incompleto)

    seed1 = rand()%0XFFFFFFFFULL;
    seed2 = rand()%0XFFFFFFFFULL;
    seed3 = rand()%0XFFFFFFFFULL;
    seed4 = rand()%0XFFFFFFFFULL;

    popsize =x->pop_size;
    chromlength = x->n_var * x->var_size;
    nintegers = chromlength / 64;
    nmask = chromlength % 64;
    if (nmask) nintegers++;

    x->chromlen = chromlength;

    poplen = x->pop_size * nintegers;

    // inizializzo i semi per la generazione dei numeri casuali

    uint64_t init[4]={seed1, seed2, seed3, seed4}, length=4;
    init_by_array64(init, length);

    for (i=nintegers; i<poplen; i++) {
    atom_setlong(&x->population[i], genrand64_int64());
    }

    if (nmask) {
        mask = UINT64_MAX << (uint64_t)(64 - nmask);
        for (i=2; i<=popsize; i++) {
            atom_setlong(&x->population[i * nintegers - 1], atom_getlong(&x->population[i
            * nintegers - 1]) & mask);
        }
    }

    BGA_gen2phen(x);
}


// ==============================================
// reinitnew: rigenero totalmente la popolazione
// ==============================================

void BGA_reinitnew(t_BGA *x)
{
    uint64_t seed1, seed2, seed3, seed4; // semi per la generazione dei numeri random
    short i; // contatore per i cicli
    short chromlength; // dimensione di un individuo (in bit)
    short nintegers; // numero di interi per rappresentare un individuo
    short poplen; // lunghezza della popolazione
    short nmask; // quantità di 1 nella mask (da sinistra)
    short popsize;
    uint64_t mask; // bitmask per la rappresentazione dell’ultimo intero (quello incompleto)

    seed1 = rand()%0XFFFFFFFFULL;
    seed2 = rand()%0XFFFFFFFFULL;
    seed3 = rand()%0XFFFFFFFFULL;
    seed4 = rand()%0XFFFFFFFFULL;

    popsize =x->pop_size;
    chromlength = x->n_var * x->var_size;
    nintegers = chromlength / 64;
    nmask = chromlength % 64;
    if (nmask) nintegers++;

    x->chromlen = chromlength;

    poplen = x->pop_size * nintegers;

    // inizializzo i semi per la generazione dei numeri casuali

    uint64_t init[4]={seed1, seed2, seed3, seed4}, length=4;
    init_by_array64(init, length);

    for (i=0; i<poplen; i++) {
        atom_setlong(&x->population[i], genrand64_int64());
    }

    if (nmask) {
        mask = UINT64_MAX << (uint64_t)(64 - nmask);
        for (i=1; i<=popsize; i++) {
            atom_setlong(&x->population[i * nintegers - 1], atom_getlong(&x->population[i
            * nintegers - 1]) & mask);
        }
    }

    BGA_gen2phen(x);
}


// =======================================================
// gen2phen: decodifico gli interi in liste di variabili
// =======================================================


void BGA_gen2phen(t_BGA *x)
{
short i, j, k, l, sum = 0;
short numvars = x->n_var;
short varsize = x->var_size;
short popsize = x->pop_size;
short nints = x->n_integers;
short count = 0;
uint64_t mask;

for (i=1; i<=popsize; i++) {
count = 0;
mask = (uint64_t)1 << 63;
for (j=0; j<numvars; j++) {
for (k = (varsize - 1); k >= 0; k--) {
sum += pow(2, k) * ((atom_getlong(&x->population[i*nints-nints+(count/64)])
& mask)? 1:0);
count++;
if (mask == 1) {
mask <<= 63;
} else mask >>= 1;
}
atom_setlong(&x->phenotype[j], sum);
sum = 0;
}
outlet_list(x->b_lambdaout, 0L, numvars, x->phenotype);
x->individual[i-1].singlefit = x->lambda_in;

l=1;
for (l=0; l<nints; l++) {
atom_setlong(&x->individual[i-1].value[l], atom_getlong(&x->population[i *
nints - nints + l]));
}
}
BGA_bubblesort(x);
}


// ======================================================
// bubblesort: ordina gli individui in base alla fitness
// (dal migliore al peggiore)
// ======================================================

void BGA_bubblesort(t_BGA *x)
{
t_cost temp;
short j,i;
short popsize = x->pop_size;

for (j=1; j < popsize; j++) {
for (i=0; i < (popsize - 1); i++) {
if (x->individual[i].singlefit < x->individual[i+1].singlefit) {
temp = x->individual[i];
x->individual[i] = x->individual[i+1];
x->individual[i+1] = temp;
}
}
}
BGA_couples(x);
}


// ================================
// couples: crea gli accoppiamenti
// ================================

void BGA_couples(t_BGA *x)
{
short i;
short m,f;
short cpopsize = x->pop_size / 2;

switch (x->matingtype) {
case 0:
for (i=1; i <= cpopsize/2; i++) {
m = (int)rand()%cpopsize;
f = (int)rand()%cpopsize;
if (m==f) {
i--;
} else {
switch (x->crosstype) {
case 0:
BGA_crossoverSP(x, m, f, i);
break;
case 1:
BGA_crossoverDP(x, m, f, i);
break;
case 2:
BGA_crossoverU(x, m, f, i);
break;
case 3:
BGA_crossoverAN(x, m, f, i);
break;
}
}
}
break;
case 1:
for (i=1; i <= cpopsize/2; i++) {
m = i * 2 - 2;
f = i * 2 - 1;
BGA_crossoverSP(x, m, f, i);
}

break;
default:
break;
}
BGA_mutation(x);
}


// ================================================
// crossoverSP: effettua il single point crossover
// ================================================

void BGA_crossoverSP(t_BGA *x, short mother, short father, short counter)
{
short i;
short nints = x->n_integers;
short chromlength = x->chromlen;
short cpopsize = x->pop_size / 2;
short splitpoint = splitpoint = (rand() % (chromlength - 1)) + 1;
short intcrossed = splitpoint / 64;
short splitinside = splitpoint % 64;
uint64_t mask1, mask2;

mask1 = (((uint64_t)1 << 63) >> splitinside) - (uint64_t)1;
mask2 = ~mask1;

for (i=0; i<nints; i++) {
if (i>intcrossed) {

atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], atom_getlong
(&x->individual[father].value[i]));
atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], atom_getlong
(&x->individual[mother].value[i]));

} else if (i == intcrossed) {

atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], (
atom_getlong(&x->individual[mother].value[i]) & mask1) | (atom_getlong(&x
->individual[father].value[i]) & mask2));
atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], (
atom_getlong(&x->individual[mother].value[i]) & mask2) | (atom_getlong(&x
->individual[father].value[i]) & mask1));

} else {

atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], atom_getlong
(&x->individual[mother].value[i]));
atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], atom_getlong
(&x->individual[father].value[i]));
}
}
}


// ================================================
// crossoverDP: effettua il double point crossover
// ================================================

void BGA_crossoverDP(t_BGA *x, short mother, short father, short counter)
{
short i,j,k;
short chromlen = x->chromlen;
short cpopsize = x->pop_size / 2;
short len = rand() % (chromlen / 2);
short initpoint = rand() % (chromlen - len);
uint64_t sum = 0;
short nints = x->n_integers;

for (i=initpoint; i<(initpoint + len); i++) {
sum += (uint64_t)(pow(2, 63 - (i%64)));
if ((i % 64) == 63) {
x->mask1[i/64] = sum;
sum = 0;
}
}

i--;
x->mask1[i/64] = sum;
sum = 0;

if ((i/64)<(nints-1)) {
for (j=(i/64)+1; j<nints; j++) {
x->mask1[j] = 0;
}
}

for (i=0; i<nints; i++) {
x->mask2[i] = ~(x->mask1[i]);
}

for (k=0; k<nints; k++) {
atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[k], (atom_getlong(&
x->individual[mother].value[k]) & x->mask1[k]) | (atom_getlong(&x->individual
[father].value[k]) & x->mask2[k]));
atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[k], (atom_getlong(&
x->individual[mother].value[k]) & x->mask2[k]) | (atom_getlong(&x->individual
[father].value[k]) & x->mask1[k]));
}
}

// ==========================================
// crossoverU: effettua il Crossover Uniform
// ==========================================

void BGA_crossoverU(t_BGA *x, short mother, short father, short counter)
{
uint64_t seed1, seed2, seed3, seed4; // semi per la generazione dei numeri random
short i;
short cpopsize = x->pop_size / 2;
short countf, countm;

short nints = x->n_integers;
uint64_t umask1[10];
uint64_t umask2[10];

countf = counter * 2 - 2 + cpopsize;
countm = counter * 2 - 1 + cpopsize;

seed1 = rand()%0XFFFFFFFFULL;
seed2 = rand()%0XFFFFFFFFULL;
seed3 = rand()%0XFFFFFFFFULL;
seed4 = rand()%0XFFFFFFFFULL;

uint64_t init[4]={seed1, seed2, seed3, seed4}, length=4;
init_by_array64(init, length);

for (i=0; i<nints; i++) {
umask1[i] = genrand64_int64();
umask2[i] = ~umask1[i];

atom_setlong(&x->individual[countf].value[i], (atom_getlong(&x->individual[mother
].value[i]) & umask1[i]) | (atom_getlong(&x->individual[father].value[i]) &
umask2[i]));
atom_setlong(&x->individual[countm].value[i], (atom_getlong(&x->individual[mother
].value[i]) & umask2[i]) | (atom_getlong(&x->individual[father].value[i]) &
umask1[i]));
}
}

// ===========================================
// crossoverAN: effettua il Crossover AND/NOT
// ===========================================

void BGA_crossoverAN(t_BGA *x, short mother, short father, short counter)
{
short i;
short cpopsize = x->pop_size / 2;
short countf, countm;

short nints = x->n_integers;

countf = counter * 2 - 2 + cpopsize;
countm = counter * 2 - 1 + cpopsize;

for (i=0; i<nints; i++) {

atom_setlong(&x->individual[countf].value[i], atom_getlong(&x->individual[mother].
value[i]) & atom_getlong(&x->individual[father].value[i]));
atom_setlong(&x->individual[countm].value[i], atom_getlong(&x->individual[mother].
value[i]) | atom_getlong(&x->individual[father].value[i]));
}
}


// ==============================================================================
// mutation: decide se avverrà (e in quale individuo) la mutazione
// e alla fine fa uscire dai vari outlet i risultati per la popolazione corrente
// ==============================================================================

void BGA_mutation(t_BGA *x)
{
short i;
short popsize = x->pop_size;
double probability = x->prob;

for (i=popsize/2; i<popsize; i++) {
if (rand()%10000 <= (probability * 100.0)) {
BGA_applyMut(x, i);
//outlet_int(x->b_out, (i*x->n_integers - x->n_integers));
//outlet_int(x->b_out, x->n_integers);
}
}

// outlets: 1. l’intera popolazione
// 2. miglior individuo
// 3. miglior fitness
// 4. attuale numero della generazione

BGA_repack(x);
}


// ========================================================================
// applyMut: effettua la mutazione sugli individui selezionati da mutation
// ========================================================================

void BGA_applyMut(t_BGA *x, short index)
{
short i;
short chromlen = x->chromlen;
short nbit = x->nbit;
long bitnumber, intnumber, tempbit = 0;
uint64_t position, bitmask = 0;
uint64_t temp;


for (i=0; i<nbit; i++) {

bitnumber = (rand()%chromlen) + 1;

if (bitnumber == tempbit) {
i--;
} else {
position = bitnumber % (uint64_t)64;
intnumber = bitnumber / 64;
bitmask = ((uint64_t)1 << (uint64_t)63) >> (position - (uint64_t)1);
temp = atom_getlong(&x->individual[index].value[intnumber]) ^ bitmask;
atom_setlong(&x->individual[index].value[intnumber], temp);
tempbit = bitnumber;
}
}
}


// =======================================================================
// repack: preleva gli interi dalla struttura (t_cost) e ricrea la lista
// degli interi con l’intera popolazione
// =======================================================================

void BGA_repack(t_BGA *x)
{
short j,k;
short popsize = x->pop_size;
short nints = x->n_integers;

for (j=1; j<=popsize; j++) {
for (k=0; k<nints; k++) {
atom_setlong(&x->population[j * nints - nints + k], atom_getlong(&x->
individual[j-1].value[k]));
}
}

outlet_int(x->b_out0, x->generationcounter);
outlet_float(x->b_out1, x->individual[0].singlefit);
outlet_list(x->b_out2, 0L, x->n_integers, x->individual[0].value);
outlet_list(x->b_out3, 0L, x->pop_size * x->n_integers, x->population);
x->generationcounter++;
}


// ========================================================
// bang: ricomincia il ciclo selezione-mutazione-crossover
// con la nuova popolazione
// ========================================================

void BGA_bang(t_BGA *x)
{
long gencount = x->generationcounter;
long refreshn = x->refreshnum;
long totref = x->totrefresh;

if (gencount % totref) {
if (gencount % refreshn) {
BGA_gen2phen(x);
} else BGA_reinit(x);
} else BGA_reinitnew(x);
}


void BGA_free(t_BGA *x)
{
if (x) {
sysmem_freeptr(x->population);
sysmem_freeptr(x->phenotype);
sysmem_freeptr(x->individual);
sysmem_freeptr(x->mask1);
sysmem_freeptr(x->mask2);
}
}


void BGA_assist(t_BGA *x, void *b, long io, long index, char *s)
{
switch (io) {

case 1:
switch (index) {
case 0:
strncpy_zero(s, "bang", 512);
break;
case 1:
strncpy_zero(s, "population size", 512);
break;
case 2:
strncpy_zero(s, "vars number", 512);
break;
case 3:
strncpy_zero(s, "var size", 512);
break;
default:
strncpy_zero(s, "lambda input", 512);
break;
}
break;

case 2:
switch (index) {
case 0:
strncpy_zero(s, "population list", 512);
break;
case 1:
strncpy_zero(s, "current best individual", 512);
break;
case 2:
strncpy_zero(s, "current best fitness", 512);
break;
case 3:
strncpy_zero(s, "generation number", 512);
break;
default:
strncpy_zero(s, "lambda output", 512);
break;
}
break;
}
}


void *BGA_new(t_symbol *s, long argc, t_atom *argv)
{
srand(time(NULL));
t_BGA *x = NULL;

x = (t_BGA *)object_alloc(BGA_class);

x->var_size = 3;
x->n_var = 8;
x->pop_size = 12;

x->population = (t_atom *)sysmem_newptr(sizeof(t_atom) * 16000);
x->phenotype = (t_atom *)sysmem_newptr(sizeof(t_atom) * 1000);
x->individual = (t_cost *)sysmem_newptr(sizeof(t_cost) * 16000);
x->mask1 = (uint64_t *)sysmem_newptr(sizeof(uint64_t) * 16);
x->mask2 = (uint64_t *)sysmem_newptr(sizeof(uint64_t) * 16);

x->matingtype = 0;
x->crosstype = 0;
x->prob = 20.0;
x->nbit = 1;
x->refreshnum = 2000;
x->totrefresh = 2000;

floatin(x, 1);

x->b_lambdaout = listout(x);
x->b_out0 = intout(x);
x->b_out1 = floatout(x);
x->b_out2 = intout(x);
x->b_out3 = listout(x);

attr_args_process(x, argc, argv);

return (x);
}
