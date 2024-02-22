#include "ext.h"
#include "ext_obex.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mt19937ar.h"

// structs

typedef struct _cost {
    double singlefit;
    t_atom* value;
} t_cost;

typedef struct _IGA

{
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
} t_IGA;

void* IGA_new(t_symbol* s, long argc, t_atom* argv);
void IGA_free(t_IGA* x);
void IGA_assist(t_IGA* x, void* b, long io, long index, char* s);
void IGA_bang(t_IGA* x);
void IGA_ft1(t_IGA* x, double fit);
void IGA_init(t_IGA* x);
void IGA_reinit(t_IGA* x);
void IGA_reinitnew(t_IGA* x);
void fillpopulation(t_IGA* x, long size);
void fillindividuals(t_IGA* x, long popsize, long nvar);
void IGA_applyfit(t_IGA* x);
void IGA_gen2phen(t_IGA* x);
void IGA_quicksort(t_IGA* x);
void IGA_couples(t_IGA* x);
void IGA_ctype(t_IGA* x, long t);
void IGA_refresh(t_IGA* x, long r);
void IGA_totrefresh(t_IGA* x, long r);
void IGA_crossoverSP(t_IGA* x, short mother, short father, short counter);
void IGA_crossoverDP(t_IGA* x, short mother, short father, short counter);
void IGA_crossoverU(t_IGA* x, short mother, short father, short counter);
void IGA_crossoverS(t_IGA* x, short mother, short father, short counter);
void IGA_crossoverA(t_IGA* x, short mother, short father, short counter);
void IGA_mutation(t_IGA* x);
void IGA_rvar(t_IGA* x, long t);
void IGA_prob(t_IGA* x, double t);
void IGA_applyMut(t_IGA* x, short index);
static int compare(const void* val1, const void* val2);

// global class pointer variable

t_class* IGA_class;

//int C74_EXPORT main(void)
void ext_main(void *r)
{
    t_class* c;

    c = class_new("bios.IGA", (method)IGA_new, (method)IGA_free, (long)sizeof(t_IGA),
        0L, A_GIMME, 0);

    class_addmethod(c, (method)IGA_assist, "assist", A_CANT, 0);
    class_addmethod(c, (method)IGA_bang, "bang", 0);
    class_addmethod(c, (method)IGA_ft1, "ft1", A_FLOAT, 0);
    class_addmethod(c, (method)IGA_init, "init", 0);
    class_addmethod(c, (method)IGA_rvar, "rvar", A_LONG, 0);
    class_addmethod(c, (method)IGA_prob, "prob", A_FLOAT, 0);
    class_addmethod(c, (method)IGA_ctype, "ctype", A_LONG, 0);
    class_addmethod(c, (method)IGA_refresh, "refresh", A_LONG, 0);
    class_addmethod(c, (method)IGA_totrefresh, "totrefresh", A_LONG, 0);

    CLASS_ATTR_LONG(c, "popsize", 0, t_IGA, pop_size);
    CLASS_ATTR_LONG(c, "numvar", 0, t_IGA, n_var);
    CLASS_ATTR_LONG_ARRAY(c, "bounds", 0, t_IGA, var_bounds, 2);

    class_register(CLASS_BOX, c); /* CLASS_NOBOX */
    IGA_class = c;

    post("∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗");
    post("−−− bios.IGA");
    post("GA Context − classical Integer GA");
    post("c 2014 July − Francesco Bianchi");
    post("v0.4.6 alpha");
    post("∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗");

    return 0;
}

/* inlets: 1. bang, init, e metodi vari
2. pop size
3. num vars
4. var size
5. lambda inlet */

void IGA_ft1(t_IGA* x, double fit)
{
    x->lambda_in = fit;
}

// ================================================
// passo il numero di bit da mutare via messaggio
// ================================================

void IGA_rvar(t_IGA* x, long t)
{
    x->rvar = t;
}

// ================================================
// passo la probabilità di mutazione via messaggio
// ================================================

void IGA_prob(t_IGA* x, double t)
{
    x->prob = t;
}

// ================================================
// passo il numero di bit da mutare via messaggio
// ================================================

void IGA_ctype(t_IGA* x, long t)
{
    x->crosstype = t;
}

// ================================================
// numero di generazioni dopo cui fare il refresh
// ================================================

void IGA_refresh(t_IGA* x, long r)
{
    x->refreshnum = r;
}

// ======================================================
// numero di generazioni dopo cui fare il refresh totale
// ======================================================

void IGA_totrefresh(t_IGA* x, long r)
{
    x->totrefresh = r;
}

// =====================================================================
// libero, ricreo e ridimensiono l’array delle coppie individui/fitness
// =====================================================================

void fillindividuals(t_IGA* x, long popsize, long nvar)
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

// ==============================================
// init: genero la popolazione iniziale (random)
// ==============================================

void IGA_init(t_IGA* x)
{
    short i, j; // contatore per i cicli
    short popsize = x->pop_size;
    short nvar = x->n_var;

    x->low = x->var_bounds[0];
    x->high = x->var_bounds[1];

    long min = x->low;
    long max = x->high;
    long modulor = (max - min) + 1;

    x->generationcounter = 0;

    if (x->rvar > x->n_var) {
        object_warn((t_object*)x, "rvar must be <= nvar\nrvar now is = nvar");
        x->rvar = x->n_var;
    }

    fillindividuals(x, x->pop_size, x->n_var);

    init_genrand(rand() % 0XFFFFFFFFUL);

    for (i = 0; i < popsize; i++) {
        for (j = 0; j < nvar; j++) {
            atom_setlong(&x->individual[i].value[j], (genrand_int32() % modulor) + min);
        }
    }
    IGA_applyfit(x);
}

// ===================================================================
// reinit: ri-genero la popolazione iniziale (random) dopo il refresh
// ===================================================================

void IGA_reinit(t_IGA* x)
{
    short i, j;
    short popsize = x->pop_size;
    short nvar = x->n_var;

    long min = x->low;
    long max = x->high;
    long modulor = (max - min) + 1;

    init_genrand(rand() % 0XFFFFFFFFUL);

    for (i = 1; i < popsize; i++) {
        for (j = 0; j < nvar; j++) {
            x->individual[i].value[j].a_w.w_long = (genrand_int32() % modulor) + min;
        }
    }
    IGA_applyfit(x);
}

// =======================================================
// reinitnew: ri-genero la popolazione iniziale (random)
// dopo il refresh totale
// =======================================================

void IGA_reinitnew(t_IGA* x)
{
    short i, j;
    short popsize = x->pop_size;
    short nvar = x->n_var;

    long min = x->low;
    long max = x->high;
    long modulor = (max - min) + 1;

    init_genrand(rand() % 0XFFFFFFFFUL);

    for (i = 0; i < popsize; i++) {
        for (j = 0; j < nvar; j++) {
            x->individual[i].value[j].a_w.w_long = (genrand_int32() % modulor) + min;
        }
    }
    IGA_applyfit(x);
}

void IGA_applyfit(t_IGA* x)
{
    short i;
    short popsize = x->pop_size;
    short nvar = x->n_var;

    for (i = 0; i < popsize; i++) {
        outlet_list(x->i_lambdaout, 0L, nvar, x->individual[i].value);
        x->individual[i].singlefit = x->lambda_in;
    }

    IGA_quicksort(x);
    IGA_couples(x);
    IGA_mutation(x);
}

// =====================================================
// quicksort: ordina gli individui in base alla fitness
// (dal migliore al peggiore)
// =====================================================

static int compare(const void* val1, const void* val2)
{
    return (((t_cost*)val1)->singlefit <= ((t_cost*)val2)->singlefit ? 1 : -1);
}

void IGA_quicksort(t_IGA* x)
{
    int length = x->pop_size;
    qsort(x->individual, length, sizeof(*x->individual), compare);
}

// ================================
// couples: crea gli accoppiamenti
// ================================

void IGA_couples(t_IGA* x)
{
    short i;
    short m, f;
    short cpopsize = x->pop_size / 2;

    switch (x->matingtype) {
    case 0:
        for (i = 1; i <= cpopsize / 2; i++) {
            m = (int)rand() % cpopsize;
            f = (int)rand() % cpopsize;
            if (m == f) {
                i--;
            } else {
                switch (x->crosstype) {
                case 0:
                    IGA_crossoverSP(x, m, f, i);
                    break;
                case 1:
                    IGA_crossoverDP(x, m, f, i);
                    break;
                case 2:
                    IGA_crossoverU(x, m, f, i);
                    break;
                case 3:
                    IGA_crossoverS(x, m, f, i);
                    break;
                case 4:
                    IGA_crossoverA(x, m, f, i);
                    break;
                }
            }
        }
        break;
    case 1:
        for (i = 1; i <= cpopsize / 2; i++) {
            m = i * 2 - 2;
            f = i * 2 - 1;
            IGA_crossoverSP(x, m, f, i);
        }
        break;
    default:
        break;
    }
}

// ==================================================
// crossoverSP: effettuare il single point crossover
// ==================================================

void IGA_crossoverSP(t_IGA* x, short mother, short father, short counter)
{
    short cpopsize = x->pop_size / 2;
    short nvar = x->n_var;
    short splitpoint = (rand() % (nvar - 1)) + 1;
    short off1 = counter * 2 - 2 + cpopsize;
    short off2 = counter * 2 - 1 + cpopsize;

    sysmem_copyptr(x->individual[mother].value, x->individual[off1].value, sizeof(t_atom*) * splitpoint);
    sysmem_copyptr(x->individual[father].value + splitpoint, x->individual[off1].value + splitpoint, sizeof(t_atom*) * (nvar - splitpoint));

    sysmem_copyptr(x->individual[father].value, x->individual[off2].value, sizeof(t_atom*) * splitpoint);
    sysmem_copyptr(x->individual[mother].value + splitpoint, x->individual[off2].value + splitpoint, sizeof(t_atom*) * (nvar - splitpoint));
}

// ==================================================
// crossoverDP: effettuare il double point crossover
// ==================================================

void IGA_crossoverDP(t_IGA* x, short mother, short father, short counter)
{
    short cpopsize = x->pop_size / 2;
    short nvar = x->n_var;
    short len = genrand_int32() % (nvar / 2);
    short splitpoint = (genrand_int32() % (nvar - len));
    short off1 = counter * 2 - 2 + cpopsize;
    short off2 = counter * 2 - 1 + cpopsize;
    short thirdpart = splitpoint + len;
    short thirddim = nvar - thirdpart;

    sysmem_copyptr(x->individual[mother].value, x->individual[off1].value, sizeof(t_atom*) * splitpoint);
    sysmem_copyptr(x->individual[father].value + splitpoint, x->individual[off1].value + splitpoint, sizeof(t_atom*) * len);

    sysmem_copyptr(x->individual[father].value, x->individual[off2].value, sizeof(t_atom*) * splitpoint);
    sysmem_copyptr(x->individual[mother].value + splitpoint, x->individual[off2].value + splitpoint, sizeof(t_atom*) * len);

    if (thirdpart < nvar) {
        sysmem_copyptr(x->individual[mother].value + thirdpart, x->individual[off1].value + thirdpart, sizeof(t_atom*) * thirddim);
        sysmem_copyptr(x->individual[father].value + thirdpart, x->individual[off2].value + thirdpart, sizeof(t_atom*) * thirddim);
    }
}

// ============================================
// crossoverU: effettuare l’ uniform crossover
// ============================================

void IGA_crossoverU(t_IGA* x, short mother, short father, short counter)
{
    short i;
    short mash;
    short cpopsize = x->pop_size / 2;
    short nvar = x->n_var;
    short off1 = counter * 2 + cpopsize - 2;

    for (i = 0; i < nvar; i++) {
        mash = rand() % 2;

        x->individual[off1 + mash].value[i].a_w.w_long = x->individual[mother].value[i].a_w.w_long;
        x->individual[off1 + (1 - mash)].value[i].a_w.w_long = x->individual[father].value[i].a_w.w_long;
        x->individual[off1 + mash].value[i].a_w.w_long = x->individual[father].value[i].a_w.w_long;
        x->individual[off1 + (1 - mash)].value[i].a_w.w_long = x->individual[mother].value[i].a_w.w_long;
    }
}

// ============================================
// crossoverS: effettuare il simple crossover
// ============================================

void IGA_crossoverS(t_IGA* x, short mother, short father, short counter)
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

void IGA_crossoverA(t_IGA* x, short mother, short father, short counter)
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

// ==============================================================================
// mutation: decide se avverrà (e in quale individuo) la mutazione
// e alla fine fa uscire dai vari outlet i risultati per la popolazione corrente
// ==============================================================================

void IGA_mutation(t_IGA* x)
{
    short i;
    short popsize = x->pop_size;
    double probability = x->prob * 100.0;

    for (i = popsize / 2; i < popsize; i++) {
        if (genrand_int32() % 10000 <= probability) {
            IGA_applyMut(x, i);
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

void IGA_applyMut(t_IGA* x, short index)
{
    short i;
    long varnumber, tempbit = 0;
    short nvar = x->n_var;
    long modulor = x->high - x->low + 1;
    long min = x->low;

    for (i = 0; i < x->rvar; i++) {

        varnumber = rand() % nvar;

        if (varnumber == tempbit) {
            i--;
        } else {
            x->individual[index].value[varnumber].a_w.w_long = (rand() % modulor) + min;
            tempbit = varnumber;
        }
    }
}

// ========================================================
// bang: ricomincia il ciclo selezione-mutazione-crossover
// con la nuova popolazione
// ========================================================

void IGA_bang(t_IGA* x)
{
    long gencount = x->generationcounter;
    long refreshn = x->refreshnum;
    long totref = x->totrefresh;

    if (gencount % totref) {
        if (gencount % refreshn) {
            IGA_applyfit(x);
        } else
            IGA_reinit(x);
    } else
        IGA_reinitnew(x);
}

void IGA_free(t_IGA* x)
{
    short i;

    if (x->individual) {
        for (i = 0; i < x->ghostpop; i++) {
            sysmem_freeptr(x->individual[i].value);
        }
        sysmem_freeptr(x->individual);
    }
}

void IGA_assist(t_IGA* x, void* b, long io, long index, char* s)
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

void* IGA_new(t_symbol* s, long argc, t_atom* argv)
{
    srand(time(NULL));
    t_IGA* x = NULL;

    x = (t_IGA*)object_alloc(IGA_class);

    x->pop_size = 12;
    x->n_var = 8;
    x->var_bounds[0] = 0;
    x->var_bounds[1] = 15;

    x->individual = (t_cost*)sysmem_newptr(sizeof(t_cost));
    x->individual[0].value = (t_atom*)sysmem_newptr(sizeof(t_atom));

    x->matingtype = 0;
    x->prob = 20.0;
    x->rvar = 1;
    x->ghostpop = 1;
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
