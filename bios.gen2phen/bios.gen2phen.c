#include "ext.h"
#include "ext_obex.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// =======================
// struttura dell’oggetto
// =======================

typedef struct _gen2phen {
    t_object ob;
    long n_var; // numero di variabili per individuo
    long var_size; // dimensione di ogni variabile (in bit)
    long n_integers; // numero di interi(int64) necessari per rappresentare un individuo
    long total_length; // lunghezza totale di un individuo (in bit)
    t_atom* individual; // array con le variabili espresse in bit
    t_atom* phenotype; // array con tutte le variabili decodificate (fenotipo)
    void* g_out; // outlet destro
    void* g_out2; // outlet sinistro
} t_gen2phen;

// ======================
// prototipi di funzione
// ======================

void* gen2phen_new(t_symbol* s, long argc, t_atom* argv);
void gen2phen_free(t_gen2phen* x);
void gen2phen_bang(t_gen2phen* x);
void gen2phen_assist(t_gen2phen* x, void* b, long io, long index, char* s);
void gen2phen_list(t_gen2phen* x, t_symbol* s, long argc, t_atom* argv);
void gen2phen_int(t_gen2phen* x, uint64_t number);

t_class* gen2phen_class; // puntatore globale alla classe

// ==============================
// inizializzazione dell’oggetto
// ==============================

int C74_EXPORT main(void)
{
    t_class* c;

    c = class_new("bios.gen2phen", (method)gen2phen_new, (method)gen2phen_free, (long)sizeof(t_gen2phen), 0L, A_GIMME, 0);

    class_addmethod(c, (method)gen2phen_bang, "bang", 0);
    class_addmethod(c, (method)gen2phen_assist, "assist", A_CANT, 0);
    class_addmethod(c, (method)gen2phen_list, "list", A_GIMME, 0);
    class_addmethod(c, (method)gen2phen_int, "int", A_LONG, 0);

    CLASS_ATTR_LONG(c, "numvar", 0, t_gen2phen, n_var);
    CLASS_ATTR_LONG(c, "varsize", 0, t_gen2phen, var_size);

    class_register(CLASS_BOX, c);
    gen2phen_class = c;

    post("∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗");
    post("−−− bios.gen2phen");
    post("GA Context − genotype to phenotype decoder");
    post("2014 July − Francesco Bianchi");
    post("v0.2.0 alpha");
    post("∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗∗");

    return 0;
}

// ===========================================
// primo inlet da sx: leggo l’intero (64 bit)
// che contiene le variabili codificate
// ===========================================

void gen2phen_int(t_gen2phen* x, uint64_t number)
{
    short i;
    uint64_t mask;
    x->n_integers = 1;
    x->total_length = x->n_var * x->var_size;

    mask = (uint64_t)1 << 63;

    for (i = 0; i < x->total_length; i++) {

        if (number & mask) {
            atom_setlong(&x->individual[i], 1);
        } else
            atom_setlong(&x->individual[i], 0);

        if (mask == 1) {
            mask = mask << 63;
        } else
            mask = mask >> 1;
    }

    gen2phen_bang(x);
}

// =================================================
// inlet sinistro: leggo la lista degli interi (64)
// che contiene le variabili codificate
// =================================================

void gen2phen_list(t_gen2phen* x, t_symbol* s, long argc, t_atom* argv)
{
    short i;
    uint64_t mask;
    x->n_integers = argc;
    x->total_length = x->n_var * x->var_size;

    mask = (uint64_t)1 << 63;

    for (i = 0; i < x->total_length; i++) {

        if (atom_getlong(&argv[i / 64]) & mask) {
            atom_setlong(&x->individual[i], 1);
        } else
            atom_setlong(&x->individual[i], 0);

        if (mask == 1) {
            mask = (uint64_t)1 << 63;
        } else
            mask = mask >> 1;
    }

    gen2phen_bang(x);
}

// ==================================================
// genero il fenotipo e lo invio all’outlet sinistro
// ==================================================

void gen2phen_bang(t_gen2phen* x)
{

    short i, j, sum = 0;
    short count = 0;
    short numvars = x->n_var;
    short varsize = x->var_size;

    for (i = 0; i < numvars; i++) {

        for (j = (varsize - 1); j >= 0; j--) {
            sum += pow(2, j) * atom_getlong(&x->individual[count]);
            count++;
        }

        atom_setlong(&x->phenotype[i], sum);
        sum = 0;
    }

    outlet_list(x->g_out2, 0L, numvars, x->phenotype);
}

void gen2phen_free(t_gen2phen* x)
{
    if (x->individual) {
        sysmem_freeptr(x->individual);
    }
    if (x->phenotype) {
        sysmem_freeptr(x->phenotype);
    }
}

void gen2phen_assist(t_gen2phen* x, void* b, long io, long index, char* s)
{
    switch (io) {
    case 1:
        strncpy_zero(s, "individual (genotype)", 512);
        break;
    case 2:
        switch (index) {
        case 0:
            strncpy_zero(s, "phenotype", 512);
            break;
        default:
            strncpy_zero(s, "outlet for future use", 512);
            break;
        }
        break;
    }
}

void* gen2phen_new(t_symbol* s, long argc, t_atom* argv)
{
    t_gen2phen* x = NULL;
    x = (t_gen2phen*)object_alloc(gen2phen_class);

    x->individual = (t_atom*)sysmem_newptr(sizeof(t_atom) * 1024);
    x->phenotype = (t_atom*)sysmem_newptr(sizeof(t_atom) * 1000);

    x->g_out = listout(x);
    x->g_out2 = listout(x);

    attr_args_process(x, argc, argv);

    return (x);
}