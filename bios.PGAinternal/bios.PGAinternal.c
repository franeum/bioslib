
#include "ext.h"
#include "ext_obex.h"
#include <math.h>
#include "mt19937ar.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXVAL 100

typedef struct _cost {
    double singlefit;
    t_atom* value;
} t_cost;

typedef struct _PGAinternal {
    t_object ob;
    long pop_size;
    long n_var;
    double lambda_in;
    long ghostpop;
    long ghostsize;
    double prob;
    long refresh;
    long totrefresh;
    long exit;
    t_cost* individual;
    long generationcounter;

    // outlets

    void* lambda_out;
    void* p_out0; // generation number
    void* p_out2; // current best individual
} t_PGAinternal;

// function prototypes

void* PGAinternal_new(t_symbol* s, long argc, t_atom* argv);
void PGAinternal_free(t_PGAinternal* x);
void PGAinternal_assist(t_PGAinternal* x, void* b, long io, long index, char* s);
//void PGAinternal_bang(t_PGAinternal* x);
//void PGAinternal_list(t_PGAinternal* x, t_symbol* s, long argc, t_atom* argv);
void PGAinternal_ft1(t_PGAinternal* x, double val);
void PGAinternal_init(t_PGAinternal* x);
void PGAinternal_applyfit(t_PGAinternal* x);
void PGAinternal_bubblesort(t_PGAinternal* x);
void PGAinternal_couples(t_PGAinternal* x);
void PGAinternal_crossover(t_PGAinternal* x, short mother, short father, short counter);
void PGAinternal_prob(t_PGAinternal* x, double t);
void PGAinternal_mutation(t_PGAinternal* x);
void PGAinternal_refresh(t_PGAinternal* x, long r);
void PGAinternal_totrefresh(t_PGAinternal* x, long r);
void PGAinternal_recreate(t_PGAinternal* x);
void PGAinternal_totrecreate(t_PGAinternal* x);
void PGAinternal_exit(t_PGAinternal* x, long r);
long fact(long n);
void fillindividual(t_PGAinternal* x, long pop, long nvar);
void rand_permutation(t_PGAinternal* x, long n, long index);

// global class pointer variable

t_class* PGAinternal_class;

void ext_main(void* r)
{
    t_class* c;

    c = class_new("bios.PGAinternal", (method)PGAinternal_new, (method)PGAinternal_free,
        (long)sizeof(t_PGAinternal), 0L, A_GIMME, 0);

    class_addmethod(c, (method)PGAinternal_assist, "assist", A_CANT, 0);
    //class_addmethod(c, (method)PGAinternal_list, "list", A_GIMME, 0);
    //class_addmethod(c, (method)PGAinternal_bang, "bang", 0);
    class_addmethod(c, (method)PGAinternal_init, "init", 0);
    class_addmethod(c, (method)PGAinternal_ft1, "ft1", A_FLOAT, 0);
    class_addmethod(c, (method)PGAinternal_prob, "prob", A_FLOAT, 0);
    class_addmethod(c, (method)PGAinternal_refresh, "refresh", A_LONG, 0);
    class_addmethod(c, (method)PGAinternal_totrefresh, "totrefresh", A_LONG, 0);
    class_addmethod(c, (method)PGAinternal_exit, "exit", A_LONG, 0);

    CLASS_ATTR_LONG(c, "popsize", 0, t_PGAinternal, pop_size);
    CLASS_ATTR_LONG(c, "numvar", 0, t_PGAinternal, n_var);

    class_register(CLASS_BOX, c); /* CLASS_NOBOX */
    PGAinternal_class = c;

    post("∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗");
    post("−−− bios.PGAinternal");
    post("GA Context − Permutations GA");
    post("c 2014 − Francesco Bianchi");
    post("v0.3.3 alpha (2014 July)");
    post("∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗");

    return 0;
}

/* inlets: 1. bang, init, e metodi vari
2. pop size
3. num vars
4. var size
5. lambda inlet */

void PGAinternal_in3(t_PGAinternal* x, long popsize)
{
    x->pop_size = popsize;
}

void PGAinternal_in2(t_PGAinternal* x, long nvar)
{
    x->n_var = nvar;
}

void PGAinternal_ft1(t_PGAinternal* x, double val)
{
    x->lambda_in = val;
}

// ================================================
// numero di generazioni dopo cui fare il refresh
// ================================================

void PGAinternal_refresh(t_PGAinternal* x, long r)
{
    x->refresh = r;
}

// ==============
// total refresh
// ==============

void PGAinternal_totrefresh(t_PGAinternal* x, long r)
{
    x->totrefresh = r;
}

void PGAinternal_exit(t_PGAinternal* x, long r)
{
    x->exit = r;
}

// ================================================
// passo la probabilità di mutazione via messaggio
// ================================================

void PGAinternal_prob(t_PGAinternal* x, double t)
{
    x->prob = t;
}

void fillindividual(t_PGAinternal* x, long pop, long nvar)
{
    short i;

    if ((pop != x->ghostpop && pop) || (nvar != x->ghostsize && nvar)) {
        if (x->individual) {
            for (i = 0; i < x->ghostpop; i++) {
                sysmem_freeptr(x->individual[i].value);
            }
            sysmem_freeptr(x->individual);

            x->individual = (t_cost*)sysmem_newptr(pop * sizeof(t_cost));

            for (i = 0; i < pop; i++) {
                x->individual[i].value = (t_atom*)sysmem_newptr(nvar * sizeof(t_atom));
            }
        } else {
            x->individual = NULL;
        }
    }
    x->ghostpop = pop;
    x->ghostsize = nvar;
}

void PGAinternal_init(t_PGAinternal* x)
{
    short i, j;
    long nvar = x->n_var;
    long popsize = x->pop_size;

    x->generationcounter = 0;

    fillindividual(x, popsize, nvar);

    for (j = 0; j < popsize; j++) {

        for (i = 0; i < nvar; i++) {
            atom_setlong(&x->individual[j].value[i], i);
        }
        rand_permutation(x, nvar, j);
    }
    PGAinternal_applyfit(x);
}

void PGAinternal_recreate(t_PGAinternal* x)
{
    short i, j;
    long nvar = x->n_var;
    long popsize = x->pop_size;

    for (j = 1; j < popsize; j++) {

        for (i = 0; i < nvar; i++) {
            atom_setlong(&x->individual[j].value[i], i);
        }
        rand_permutation(x, nvar, j);
    }
    PGAinternal_applyfit(x);
}

void PGAinternal_totrecreate(t_PGAinternal* x)
{
    short i, j;
    long nvar = x->n_var;
    long popsize = x->pop_size;

    for (j = 0; j < popsize; j++) {

        for (i = 0; i < nvar; i++) {
            atom_setlong(&x->individual[j].value[i], i);
        }
        rand_permutation(x, nvar, j);
    }
    PGAinternal_applyfit(x);
}

void PGAinternal_applyfit(t_PGAinternal* x)
{
    short i;
    long nvar = x->n_var;
    long popsize = x->pop_size;

    for (i = 0; i < popsize; i++) {
        outlet_list(x->lambda_out, 0L, nvar, x->individual[i].value);
        x->individual[i].singlefit = x->lambda_in;
    }
    PGAinternal_bubblesort(x);
}

// =================================================================================
// bubblesort: ordina gli individui in base alla fitness (dal migliore al peggiore)
// =================================================================================

void PGAinternal_bubblesort(t_PGAinternal* x)
{
    t_cost temp;
    short j, i;
    short popsize = x->pop_size;

    for (j = 1; j < popsize; j++) {
        for (i = 0; i < (popsize - 1); i++) {
            if (x->individual[i].singlefit <= x->individual[i + 1].singlefit) {
                temp = x->individual[i];
                x->individual[i] = x->individual[i + 1];
                x->individual[i + 1] = temp;
            }
        }
    }
    PGAinternal_couples(x);
}

// ========================================
// couples: crea gli accoppiamenti
// ========================================

void PGAinternal_couples(t_PGAinternal* x)
{
    short i;
    short m, f;
    short cpopsize = x->pop_size / 2;

    for (i = 1; i <= cpopsize / 2; i++) {
        m = (int)rand() % cpopsize;
        f = (int)rand() % cpopsize;
        if (m == f) {
            i--;
        } else
            PGAinternal_crossover(x, m, f, i);
    }
    PGAinternal_mutation(x);
}

// ==================================================
// crossover: effettuare il crossover
// ==================================================

void PGAinternal_crossover(t_PGAinternal* x, short mother, short father, short counter)
{
    short i, j;
    short cpopsize = x->pop_size / 2;
    short nvar = x->n_var;
    short len = rand() % (nvar / 2);
    short splitpoint = (rand() % (nvar - len));

    long temp;
    short determ;
    short index;

    for (i = splitpoint; i < len + splitpoint; i++) {
        atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], atom_getlong(&x->individual[mother].value[i]));
        atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], atom_getlong(&x->individual[father].value[i]));
    }

    index = 0;

    for (i = 0; i < nvar; i++) {

        determ = 0;
        temp = atom_getlong(&x->individual[father].value[i]);

        for (j = splitpoint; j < len + splitpoint; j++)
            if (temp == atom_getlong(&x->individual[mother].value[j]))
                determ++;

        if (determ == 0 && index < splitpoint) {
            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[index], temp);
            index++;
        } else if (determ == 0 && index >= splitpoint) {
            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[index + len],
                temp);
            index++;
        }
    }
    index = 0;

    for (i = 0; i < nvar; i++) {

        determ = 0;
        temp = atom_getlong(&x->individual[mother].value[i]);

        for (j = splitpoint; j < len + splitpoint; j++)
            if (temp == atom_getlong(&x->individual[father].value[j]))
                determ++;

        if (determ == 0 && index < splitpoint) {
            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[index], temp);
            index++;
        } else if (determ == 0 && index >= splitpoint) {
            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[index + len],
                temp);
            index++;
        }
    }
}

// ==============================================================================
// mutation: decide se avverrà (e in quale individuo) la mutazione
// e alla fine fa uscire dai vari outlet i risultati per la popolazione corrente
// ==============================================================================

void PGAinternal_mutation(t_PGAinternal* x)
{
    short i;
    short popsize = x->pop_size;
    long nvar = x->n_var;
    double probability = x->prob;
    double currentbestfitness;
    long exit = x->exit;
    long refresh = x->refresh;
    long totrefresh = x->totrefresh;

    for (i = popsize / 2; i < popsize; i++) {
        if (rand() % 10000 <= (probability * 100.0)) {
            rand_permutation(x, nvar, i);
        }
    }

    /* outlets: 1. l’intera popolazione
    2. miglior individuo
    3. miglior fitness
    4. attuale numero della generazione */

    currentbestfitness = x->individual[0].singlefit;

    if (currentbestfitness == 1.) {
        outlet_int(x->p_out0, x->generationcounter);
        outlet_list(x->p_out2, 0L, x->n_var, x->individual[0].value);
    } else if (currentbestfitness != 1. && (x->generationcounter < exit)) {
        x->generationcounter++;
        if (x->generationcounter % totrefresh) {
            if (x->generationcounter % refresh) {
                PGAinternal_applyfit(x);
            } else
                PGAinternal_recreate(x);
        } else
            PGAinternal_totrecreate(x);
    } else
        post("no solution found!! retry!");
}

void PGAinternal_free(t_PGAinternal* x)
{
    short i;

    if (x->individual) {
        for (i = 0; i < x->ghostpop; i++) {
            sysmem_freeptr(x->individual[i].value);
        }
        sysmem_freeptr(x->individual);
    }
}

void PGAinternal_assist(t_PGAinternal* x, void* b, long io, long index, char* s)
{
    switch (io) {

    case 1:
        switch (index) {
        case 0:
            strncpy_zero(s, "messages", 512);
            break;
        default:
            strncpy_zero(s, "lambda input", 512);
            break;
        }
        break;

    case 2:
        switch (index) {
        case 0:
            strncpy_zero(s, "solution", 512);
            break;
        case 1:
            strncpy_zero(s, "generation number", 512);
            break;
        default:
            strncpy_zero(s, "lambda output", 512);
            break;
        }
        break;
    }
}

void* PGAinternal_new(t_symbol* s, long argc, t_atom* argv)
{
    srand(time(NULL));
    t_PGAinternal* x = NULL;

    x = (t_PGAinternal*)object_alloc(PGAinternal_class);

    x->individual = (t_cost*)sysmem_newptr(sizeof(t_cost));
    x->individual[0].value = (t_atom*)sysmem_newptr(sizeof(t_atom));

    x->refresh = 1000;
    x->totrefresh = 2000;
    x->exit = 5000;
    x->pop_size = 12;
    x->n_var = 8;

    floatin(x, 1);
    x->lambda_out = listout(x);
    x->p_out0 = intout(x);
    x->p_out2 = listout(x);

    attr_args_process(x, argc, argv);

    return (x);
}

void rand_permutation(t_PGAinternal* x, long n, long index)
{
    if (n > 1) {
        long i;
        for (i = 0; i < n - 1; i++) {
            long j = i + rand() / (RAND_MAX / (n - i) + 1);
            long t = atom_getlong(&x->individual[index].value[j]);
            atom_setlong(&x->individual[index].value[j], atom_getlong(&x->individual[index].value[i]));
            atom_setlong(&x->individual[index].value[i], t);
        }
    }
}
