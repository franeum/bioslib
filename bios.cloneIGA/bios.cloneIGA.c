#include "ext.h"
#include "ext_obex.h"
#include <math.h>
#include "mt19937ar.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// structs

typedef struct _cost {
    double singlefit;
    t_atom* value;
} t_cost;

typedef struct _cloneIGA {
    t_object ob;
    long pop_size; // dimensione della popolazione
    long n_var; // numero di variabili per individuo
    long low; // valore minimo di una variabile
    long high; // valore massimo di una variabile
    long var_bounds[2]; // array per i limiti di una variabile
    short ghostpop;
    long totaldim; // dimensione totale della matrice (n_var * pop_size)
    double lambda_in; // valore di fitness in entrata dal lambda inlet
    short matingtype; // tipo di algoritmo per l’accoppiamento
    short crosstype;
    long fraglen;
    long changex;
    long muttype;
    long mutstepsize;
    long refreshnum; // numero di generazioni dopo le quali fare il refresh
    long totrefresh; // numero di generazioni dopo le quali fare il refresh totale
    long rvar; // numero di bit da mutare per ogni individuo
    double prob; // probabilità di mutazione
    t_cost* individual; // array con le coppie individui/fitness
    long generationcounter;

    // outlets

    void* i_lambdaout;
    void* i_out0;
    void* i_out1;
    void* i_out2;
    void* i_out3;
} t_cloneIGA;

// prototipi

void* cloneIGA_new(t_symbol* s, long argc, t_atom* argv);
void cloneIGA_free(t_cloneIGA* x);
void cloneIGA_assist(t_cloneIGA* x, void* b, long io, long index, char* s);
void cloneIGA_bang(t_cloneIGA* x);
void cloneIGA_ft1(t_cloneIGA* x, double fit);
void cloneIGA_list(t_cloneIGA* x, t_symbol* s, long argc, t_atom* argv);
void cloneIGA_init(t_cloneIGA* x);
void cloneIGA_reinit(t_cloneIGA* x);
void cloneIGA_reinitnew(t_cloneIGA* x);
void fillpopulation(t_cloneIGA* x, long size);
void fillindividuals(t_cloneIGA* x, long popsize, long nvar);
void cloneIGA_applyfit(t_cloneIGA* x);
void cloneIGA_gen2phen(t_cloneIGA* x);
void cloneIGA_bubblesort(t_cloneIGA* x);
void cloneIGA_couples(t_cloneIGA* x);
void cloneIGA_ctype(t_cloneIGA* x, long t);
void cloneIGA_refresh(t_cloneIGA* x, long r);
void cloneIGA_totrefresh(t_cloneIGA* x, long r);
void cloneIGA_fraglen(t_cloneIGA* x, long t);
void cloneIGA_changexover(t_cloneIGA* x, long t);
void cloneIGA_muttype(t_cloneIGA* x, long t);
void cloneIGA_mutstep(t_cloneIGA* x, long t);
void cloneIGA_crossoverSP(t_cloneIGA* x, short mother, short father, short counter);
void cloneIGA_crossoverDP(t_cloneIGA* x, short mother, short father, short counter);
void cloneIGA_crossoverU(t_cloneIGA* x, short mother, short father, short counter);
void cloneIGA_crossoverS(t_cloneIGA* x, short mother, short father, short counter);
void cloneIGA_crossoverA(t_cloneIGA* x, short mother, short father, short counter);
void cloneIGA_crossoverLC(t_cloneIGA* x, short mother, short father, short counter);
void cloneIGA_mutation(t_cloneIGA* x);
void cloneIGA_rvar(t_cloneIGA* x, long t);
void cloneIGA_prob(t_cloneIGA* x, double t);
void cloneIGA_applyMut(t_cloneIGA* x, short index);
void cloneIGA_applyMutStep(t_cloneIGA* x, short index);

// global class pointer variable

t_class* cloneIGA_class;

void ext_main(void* r)
{
    t_class* c;

    c = class_new("bios.cloneIGA", (method)cloneIGA_new, (method)cloneIGA_free, (long)sizeof(t_cloneIGA), 0L, A_GIMME, 0);

    class_addmethod(c, (method)cloneIGA_assist, "assist", A_CANT, 0);
    class_addmethod(c, (method)cloneIGA_bang, "bang", 0);
    class_addmethod(c, (method)cloneIGA_list, "list", A_GIMME, 0);
    class_addmethod(c, (method)cloneIGA_ft1, "ft1", A_FLOAT, 0);
    class_addmethod(c, (method)cloneIGA_init, "init", 0);
    class_addmethod(c, (method)cloneIGA_rvar, "rvar", A_LONG, 0);
    class_addmethod(c, (method)cloneIGA_prob, "prob", A_FLOAT, 0);
    class_addmethod(c, (method)cloneIGA_ctype, "ctype", A_LONG, 0);
    class_addmethod(c, (method)cloneIGA_fraglen, "fraglen", A_LONG, 0);
    class_addmethod(c, (method)cloneIGA_changexover, "changex", A_LONG, 0);
    class_addmethod(c, (method)cloneIGA_muttype, "mutx", A_LONG, 0);
    class_addmethod(c, (method)cloneIGA_mutstep, "mutstep", A_LONG, 0);
    class_addmethod(c, (method)cloneIGA_refresh, "refresh", A_LONG, 0);
    class_addmethod(c, (method)cloneIGA_totrefresh, "totrefresh", A_LONG, 0);

    CLASS_ATTR_LONG(c, "popsize", 0, t_cloneIGA, pop_size);
    CLASS_ATTR_LONG_ARRAY(c, "bounds", 0, t_cloneIGA, var_bounds, 2);

    class_register(CLASS_BOX, c); /* CLASS_NOBOX */
    cloneIGA_class = c;

    post("∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗");
    post("−−− bios.cloneIGA");
    post("GA Context − cloned Integer GA");
    post("2014 − Francesco Bianchi");
    post("v0.1.1 alpha (2014 July)");
    post("∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗");

    return 0;
}

/* inlets: 1. bang, init, e metodi vari
2. pop size
3. num vars
4. var size
5. lambda inlet */

void cloneIGA_ft1(t_cloneIGA* x, double fit)
{
    x->lambda_in = fit;
}

void cloneIGA_list(t_cloneIGA* x, t_symbol* s, long argc, t_atom* argv)
{
    x->n_var = argc;
    long popsize = x->pop_size;
    long nvar = x->n_var;
    short i, j;

    x->low = x->var_bounds[0];
    x->high = x->var_bounds[1];

    x->generationcounter = 0;

    fillindividuals(x, popsize, nvar);

    for (i = 0; i < popsize; i++) {
        for (j = 0; j < nvar; j++) {
            atom_setlong(&x->individual[i].value[j], atom_getlong(&argv[j]));
        }
    }
    cloneIGA_applyfit(x);
}

void cloneIGA_mutstep(t_cloneIGA* x, long t)
{
    x->mutstepsize = t;
}

void cloneIGA_muttype(t_cloneIGA* x, long t)
{
    x->muttype = t;
}

void cloneIGA_fraglen(t_cloneIGA* x, long t)
{
    x->fraglen = t;
}

void cloneIGA_changexover(t_cloneIGA* x, long t)
{
    x->changex = t;
}

// ================================================
// passo il numero di bit da mutare via messaggio
// ================================================

void cloneIGA_rvar(t_cloneIGA* x, long t)
{
    x->rvar = t;
}

// ================================================
// passo la probabilità di mutazione via messaggio
// ================================================

void cloneIGA_prob(t_cloneIGA* x, double t)
{
    x->prob = t;
}

// ================================================
// passo il numero di bit da mutare via messaggio
// ================================================

void cloneIGA_ctype(t_cloneIGA* x, long t)
{
    x->crosstype = t;
}

// ================================================
// numero di generazioni dopo cui fare il refresh
// ================================================

void cloneIGA_refresh(t_cloneIGA* x, long r)
{
    x->refreshnum = r;
}

// ======================================================
// numero di generazioni dopo cui fare il refresh totale
// ======================================================

void cloneIGA_totrefresh(t_cloneIGA* x, long r)
{
    x->totrefresh = r;
}

// =====================================================================
// libero, ricreo e ridimensiono l’array delle coppie individui/fitness
// =====================================================================

void fillindividuals(t_cloneIGA* x, long popsize, long nvar)
{
    short i;

    if ((((x->totaldim / popsize) != nvar) || ((x->totaldim / nvar) != popsize)) && popsize && nvar) {
        if (x->individual) {
            for (i = 0; i < x->ghostpop; i++) {
                sysmem_freeptr(x->individual[i].value);
            }
            sysmem_freeptr(x->individual);

            x->individual = (t_cost*)sysmem_newptr(popsize * sizeof(t_cost));

            for (i = 0; i < popsize; i++) {
                x->individual[i].value = (t_atom*)sysmem_newptr(nvar * sizeof(t_atom));
            }
        } else {
            x->individual = NULL;
            x->individual[0].value = NULL;
        }
    }
    x->totaldim = x->pop_size * x->n_var;
    x->ghostpop = popsize;
}

// ============================================================
// con il metodo init genero la popolazione iniziale (random)
// ============================================================

void cloneIGA_init(t_cloneIGA* x)
{
    unsigned long seed1, seed2, seed3, seed4; // semi per la generazione dei numeri random
    short i, j; // contatore per i cicli
    short popsize = x->pop_size;
    short nvar = x->n_var;

    x->low = x->var_bounds[0];
    x->high = x->var_bounds[1];

    long min = x->low;
    long max = x->high;
    long modulor = (max - min) + 1;

    x->generationcounter = 0;

    seed1 = rand() % 0XFFFFFFFFUL;
    seed2 = rand() % 0XFFFFFFFFUL;
    seed3 = rand() % 0XFFFFFFFFUL;
    seed4 = rand() % 0XFFFFFFFFUL;

    if (x->rvar > x->n_var) {
        object_warn((t_object*)x, "rvar must be <= nvar\nrvar now is = nvar");
        x->rvar = x->n_var;
    }

    // inizializzo i semi per la generazione dei numeri casuali

    fillindividuals(x, x->pop_size, x->n_var);

    unsigned long init[4] = { seed1, seed2, seed3, seed4 }, length = 4;
    init_by_array(init, length);

    for (i = 0; i < popsize; i++) {
        for (j = 0; j < nvar; j++) {
            atom_setlong(&x->individual[i].value[j], (genrand_int32() % modulor) + min);
        }
    }
    cloneIGA_applyfit(x);
}

// ==================================================================
// reinit ri-genero la popolazione iniziale (random) dopo il refresh
// ==================================================================

void cloneIGA_reinit(t_cloneIGA* x)
{
    unsigned long seed1, seed2, seed3, seed4; // semi per la generazione dei numeri random
    short i, j; // contatore per i cicli
    short popsize = x->pop_size;
    short nvar = x->n_var;

    long min = x->low;
    long max = x->high;
    long modulor = (max - min) + 1;

    seed1 = rand() % 0XFFFFFFFFUL;
    seed2 = rand() % 0XFFFFFFFFUL;
    seed3 = rand() % 0XFFFFFFFFUL;
    seed4 = rand() % 0XFFFFFFFFUL;

    unsigned long init[4] = { seed1, seed2, seed3, seed4 }, length = 4;
    init_by_array(init, length);

    for (i = 1; i < popsize; i++) {
        for (j = 0; j < nvar; j++) {
            atom_setlong(&x->individual[i].value[j], (genrand_int32() % modulor) + min);
        }
    }
    cloneIGA_applyfit(x);
}

// ======================================================
// reinitnew: ri-genero la popolazione iniziale (random)
// dopo il refresh totale
// ======================================================

void cloneIGA_reinitnew(t_cloneIGA* x)
{
    unsigned long seed1, seed2, seed3, seed4; // semi per la generazione dei numeri random
    short i, j; // contatore per i cicli
    short popsize = x->pop_size;
    short nvar = x->n_var;
    long min = x->low;
    long max = x->high;
    long modulor = (max - min) + 1;

    seed1 = rand() % 0XFFFFFFFFUL;
    seed2 = rand() % 0XFFFFFFFFUL;
    seed3 = rand() % 0XFFFFFFFFUL;
    seed4 = rand() % 0XFFFFFFFFUL;

    unsigned long init[4] = { seed1, seed2, seed3, seed4 }, length = 4;
    init_by_array(init, length);

    for (i = 0; i < popsize; i++) {
        for (j = 0; j < nvar; j++) {
            atom_setlong(&x->individual[i].value[j], (genrand_int32() % modulor) + min);
        }
    }
    cloneIGA_applyfit(x);
}

void cloneIGA_applyfit(t_cloneIGA* x)
{
    short i;
    short popsize = x->pop_size;
    short nvar = x->n_var;

    for (i = 0; i < popsize; i++) {
        outlet_list(x->i_lambdaout, 0L, nvar, x->individual[i].value);
        x->individual[i].singlefit = x->lambda_in;
    }
    cloneIGA_bubblesort(x);
    cloneIGA_couples(x);
    cloneIGA_mutation(x);
}

// ======================================================
// bubblesort: ordina gli individui in base alla fitness
// (dal migliore al peggiore)
// ======================================================

void cloneIGA_bubblesort(t_cloneIGA* x)
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
}

// ================================
// couples: crea gli accoppiamenti
// ================================

void cloneIGA_couples(t_cloneIGA* x)
{
    short i;
    short m, f;
    short cpopsize = x->pop_size / 2;

    for (i = 1; i <= cpopsize / 2; i++) {
        m = (int)rand() % cpopsize;
        f = (int)rand() % cpopsize;
        if (m == f) {
            i--;
        } else {

            if (x->generationcounter < x->changex) {
                cloneIGA_crossoverLC(x, m, f, i);
            } else {

                switch (x->crosstype) {
                case 0:
                    cloneIGA_crossoverSP(x, m, f, i);
                    break;
                case 1:
                    cloneIGA_crossoverDP(x, m, f, i);
                    break;
                case 2:
                    cloneIGA_crossoverU(x, m, f, i);
                    break;
                case 3:
                    cloneIGA_crossoverS(x, m, f, i);
                    break;
                case 4:
                    cloneIGA_crossoverA(x, m, f, i);
                    break;
                }
            }
        }
    }
}

// ==================================================
// crossoverSP: effettuare il single point crossover
// ==================================================

void cloneIGA_crossoverSP(t_cloneIGA* x, short mother, short father, short counter)
{

    short i;
    short cpopsize = x->pop_size / 2;
    short nvar = x->n_var;
    short splitpoint = (rand() % (nvar - 1)) + 1;

    for (i = 0; i < nvar; i++) {
        if (i <= splitpoint) {
            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], atom_getlong(&x->individual[mother].value[i]));
            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], atom_getlong(&x->individual[father].value[i]));
        } else {
            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], atom_getlong(&x->individual[father].value[i]));
            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], atom_getlong(&x->individual[mother].value[i]));
        }
    }
}

// ==================================================
// crossoverDP: effettuare il double point crossover
// ==================================================

void cloneIGA_crossoverDP(t_cloneIGA* x, short mother, short father, short counter)
{
    short i;
    short cpopsize = x->pop_size / 2;
    short nvar = x->n_var;
    short len = rand() % (nvar / 2);
    short splitpoint = (rand() % (nvar - len));

    for (i = 0; i < nvar; i++) {
        if (i < splitpoint || i > (splitpoint + len)) {
            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], atom_getlong(&x->individual[mother].value[i]));
            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], atom_getlong(&x->individual[father].value[i]));
        } else {
            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], atom_getlong(&x->individual[father].value[i]));
            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], atom_getlong(&x->individual[mother].value[i]));
        }
    }
}

// ============================================
// crossoverU: effettuare l’ uniform crossover
// ============================================

void cloneIGA_crossoverU(t_cloneIGA* x, short mother, short father, short counter)
{
    short i;
    short mash;
    short cpopsize = x->pop_size / 2;
    short nvar = x->n_var;

    for (i = 0; i < nvar; i++) {
        mash = rand() % 2;
        if (mash == 0) {
            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], atom_getlong(&x->individual[mother].value[i]));
            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], atom_getlong(&x->individual[father].value[i]));
        } else {
            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], atom_getlong(&x->individual[father].value[i]));
            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], atom_getlong(&x->individual[mother].value[i]));
        }
    }
}

// ============================================
// crossoverS: effettuare il simple crossover
// ============================================

void cloneIGA_crossoverS(t_cloneIGA* x, short mother, short father, short counter)
{
    short i;
    short cpopsize = x->pop_size / 2;
    short nvar = x->n_var;
    short splitpoint = (rand() % (nvar - 1)) + 1;
    double alfa = (rand() % 100) / 100.0;

    for (i = 0; i < nvar; i++) {
        if (i <= splitpoint) {
            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], atom_getlong(&x->individual[mother].value[i]));
            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], atom_getlong(&x->individual[father].value[i]));
        } else {
            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], (alfa * atom_getlong(&x->individual[father].value[i]) + ((1 - alfa) * atom_getlong(&x->individual[mother].value[i]))));
            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], (alfa * atom_getlong(&x->individual[mother].value[i]) + ((1 - alfa) * atom_getlong(&x->individual[father].value[i]))));
        }
    }
}

// ============================================
// crossoverA: effettuare l’Average crossover
// ============================================

void cloneIGA_crossoverA(t_cloneIGA* x, short mother, short father, short counter)
{
    short i;
    short cpopsize = x->pop_size / 2;
    short nvar = x->n_var;
    short splitpoint = (rand() % (nvar - 1)) + 1;

    for (i = 0; i < nvar; i++) {
        if (i <= splitpoint) {
            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], atom_getlong(&x->individual[mother].value[i]));
            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], (atom_getlong(&x->individual[father].value[i]) + atom_getlong(&x->individual[father].value[i])) / 2);
        } else {
            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], (atom_getlong(&x->individual[father].value[i]) + atom_getlong(&x->individual[mother].value[i])) / 2);
            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], sqrt((atom_getlong(&x->individual[father].value[i]) * atom_getlong(&x->individual[mother].value[i]))));
        }
    }
}

// =======================
// Locus Change Crossover
// =======================

void cloneIGA_crossoverLC(t_cloneIGA* x, short mother, short father, short counter)
{
    short i;
    long fragmentlen = x->fraglen;
    short cpopsize = x->pop_size / 2;
    short nvar = x->n_var;
    long point1, point2;
    point1 = rand() % (nvar - fragmentlen);
    point2 = rand() % (nvar - fragmentlen);

    for (i = 0; i < nvar; i++) {

        if (i < point1 || i >= (point1 + fragmentlen)) {

            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i], atom_getlong(&x->individual[mother].value[i]));
        } else
            atom_setlong(&x->individual[counter * 2 - 2 + cpopsize].value[i],
                atom_getlong(&x->individual[father].value[point2 + (i - point1)]));

        if (i < point2 || i >= (point2 + fragmentlen)) {

            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i], atom_getlong(&x->individual[father].value[i]));
        } else
            atom_setlong(&x->individual[counter * 2 - 1 + cpopsize].value[i],
                atom_getlong(&x->individual[mother].value[point1 + (i - point2)]));
    }
}

// ================================================================
// mutation: decide se avverrà (e in quale individuo) la mutazione
// ================================================================

void cloneIGA_mutation(t_cloneIGA* x)
{
    short i;
    short popsize = x->pop_size;
    double probability = x->prob;

    for (i = popsize / 2; i < popsize; i++) {
        if (rand() % 10000 <= (probability * 100.0)) {
            switch (x->muttype) {

            case 0:
                cloneIGA_applyMut(x, i);
                break;

            default:
                cloneIGA_applyMutStep(x, i);
                break;
            }
        }
    }

    /* outlets: 1. l’intera popolazione
    2. miglior individuo
    3. miglior fitness
    4. attuale numero della generazione */

    outlet_int(x->i_out0, x->generationcounter);
    outlet_float(x->i_out1, x->individual[0].singlefit);
    outlet_list(x->i_out2, 0L, x->n_var, x->individual[0].value);
    x->generationcounter++;
}

// ========================================================================
// applyMut: effettua la mutazione sugli individui selezionati da mutation
// ========================================================================

void cloneIGA_applyMut(t_cloneIGA* x, short index)
{
    short i;
    long varnumber, tempbit = 0;
    short nvar = x->n_var;

    // mutazione tradizionale

    for (i = 0; i < x->rvar; i++) {

        varnumber = rand() % nvar;

        if (varnumber == tempbit) {
            i--;
        } else {

            atom_setlong(&x->individual[index].value[varnumber], (rand() % (x->high - x->low + 1)) + x->low);
            tempbit = varnumber;
        }
    }
}

void cloneIGA_applyMutStep(t_cloneIGA* x, short index)
{
    short i;
    long varnumber, tempbit = 0;
    short nvar = x->n_var;
    long tempval;
    long min = x->low;
    long max = x->high;
    long stepsize = x->mutstepsize;
    long stepmodule = stepsize * 2 + 1;

    for (i = 0; i < x->rvar; i++) {

        varnumber = rand() % nvar;

        if (varnumber == tempbit) {
            i--;
        } else {
            tempval = ((rand() % stepmodule) - stepsize) + atom_getlong(&x->individual[index].value[varnumber]);
            if (tempval < min) {
                atom_setlong(&x->individual[index].value[varnumber], min);
            } else if (tempval > max) {
                atom_setlong(&x->individual[index].value[varnumber], min);
            } else
                atom_setlong(&x->individual[index].value[varnumber], tempval);
            tempbit = varnumber;
        }
    }
}

// ========================================================
// bang: ricomincia il ciclo selezione-mutazione-crossover
// con la nuova popolazione
// ========================================================

void cloneIGA_bang(t_cloneIGA* x)
{
    long gencount = x->generationcounter;
    long refreshn = x->refreshnum;
    long totref = x->totrefresh;

    if (gencount % totref) {
        if (gencount % refreshn) {
            cloneIGA_applyfit(x);
        } else
            cloneIGA_reinit(x);
    } else
        cloneIGA_reinitnew(x);
}

void cloneIGA_free(t_cloneIGA* x)
{
    short i;

    if (x->individual) {
        for (i = 0; i < x->ghostpop; i++) {
            sysmem_freeptr(x->individual[i].value);
        }
        sysmem_freeptr(x->individual);
    }
}

void cloneIGA_assist(t_cloneIGA* x, void* b, long io, long index, char* s)
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
            strncpy_zero(s, "min var", 512);
            break;
        case 4:
            strncpy_zero(s, "max var", 512);
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

void* cloneIGA_new(t_symbol* s, long argc, t_atom* argv)
{

    srand(time(NULL));
    t_cloneIGA* x = NULL;

    x = (t_cloneIGA*)object_alloc(cloneIGA_class);
    x->pop_size = 12;
    x->var_bounds[0] = 0;
    x->var_bounds[1] = 15;

    x->individual = (t_cost*)sysmem_newptr(sizeof(t_cost));
    x->individual[0].value = (t_atom*)sysmem_newptr(sizeof(t_atom));

    x->matingtype = 0;
    x->prob = 20.0;
    x->rvar = 1;
    x->ghostpop = 1;
    x->changex = 1;
    x->mutstepsize = 2;
    x->fraglen = 1;
    x->muttype = 0;
    x->refreshnum = 2000;
    x->totrefresh = 2000;

    floatin(x, 1);

    x->i_lambdaout = listout(x);
    x->i_out0 = intout(x);
    x->i_out1 = floatout(x);
    x->i_out2 = listout(x);
    x->i_out3 = listout(x);

    attr_args_process(x, argc, argv);

    return (x);
}
