#define ACCELERATE_NEW_LAPACK
#define ACCELERATE_LAPACK_ILP64
#include "ext.h"			// standard Max include, always required (except in Jitter)
#include "ext_obex.h"		// required for "new" style objects
#include "z_dsp.h"			// required for MSP objects
#include<simd/simd.h>
#include<Accelerate/Accelerate.h>

typedef SparseMatrix_Double t_chebyshev;

C74_HIDDEN void chebyshev_init(t_chebyshev*const matrix, long const count) {
    long const nz = ( count / 2 + 1 ) * ( ( count + 1 ) / 2 );
    void * const memory = sysmem_resizeptrclear(matrix->data, ( count + 1 ) * sizeof(long) + nz * (sizeof(double) + sizeof(int)));
    matrix->structure.rowCount = (int const)count;
    matrix->structure.columnCount = (int const)count;
    matrix->structure.rowIndices = (int*const)(memory + nz * (sizeof(double)));
    matrix->structure.columnStarts = (long*const)(memory + nz * (sizeof(double) + sizeof(int)));
    matrix->structure.attributes = (SparseAttributes_t const) {
        .transpose = 0,
        .triangle = SparseUpperTriangle,
        .kind = SparseTriangular,
        ._reserved = 0,
        ._allocatedBySparse = 0,
    };
    matrix->structure.blockSize = 1;
    matrix->data = (double*const)memory;
    if ( count < 2 ) {
        matrix->data[0] = 1;
        matrix->structure.rowIndices[0] = 0;
        matrix->structure.columnStarts[0] = 0;
        matrix->structure.columnStarts[1] = 1;
    } else {
        vDSP_vclrD(matrix->data, 1, nz);
        matrix->data[0] = 1;
        matrix->data[1] = 2;
        matrix->structure.rowIndices[0] = 0;
        matrix->structure.rowIndices[1] = 1;
        matrix->structure.columnStarts[0] = 0;
        matrix->structure.columnStarts[1] = 1;
        matrix->structure.columnStarts[2] = 2;
        for ( unsigned long c = 2, C = count ; c < C ; ++ c ) {
            unsigned long const j = matrix->structure.columnStarts[c-2];
            unsigned long const k = matrix->structure.columnStarts[c-1];
            unsigned long const l = matrix->structure.columnStarts[c-0];
            unsigned long const p = c % 2;
            unsigned long const q = c / 2 + 1;
            for ( register unsigned long r = 0, R = q ; r < R ; ++ r )
                matrix->structure.rowIndices[l+r] = (int const)(2 * r + p);
            vDSP_vsmulD(matrix->data + k, 1, (double const[]){2.0}, matrix->data + l + 1 - p, 1, l - k);
            vDSP_vsubD(matrix->data + j, 1, matrix->data + l, 1, matrix->data + l, 1, k - j);
            matrix->structure.columnStarts[c+1] = l + q;
        }
    }
}
C74_HIDDEN void chebyshev_free(t_chebyshev*const matrix) {
    sysmem_freeptr(matrix->data);
    matrix->data = NULL;
}
// Col-Vector <- Col-Vector
C74_HIDDEN void chebyshev_mv(t_chebyshev const*const matrix, double*const y, double*const x) {
    SparseMultiply(*matrix,
                   (DenseVector_Double const) { .data = x, .count = matrix->structure.columnCount, },
                   (DenseVector_Double const) { .data = y, .count = matrix->structure.rowCount, });
}
// Row-Major <- Sp * Row-Major
C74_HIDDEN void chebyshev_mm(t_chebyshev const*const matrix, double*const y, long const ldy, double*const x,long const ldx, long const length) {
    SparseMultiply(*matrix,
                   (DenseMatrix_Double const) {
        .rowCount = (int const)length,
        .columnCount = matrix->structure.columnCount,
        .columnStride = (int const)ldx,
        .attributes = (SparseAttributes_t const) {
            .transpose = 1,
            .triangle = 0,
            .kind = SparseOrdinary,
            ._reserved = 0,
            ._allocatedBySparse = 0,
        },
        .data = x,
    },
                   (DenseMatrix_Double const) {
        .rowCount = (int const)length,
        .columnCount = matrix->structure.rowCount,
        .columnStride = (int const)ldy,
        .attributes = (SparseAttributes_t const) {
            .transpose = 1,
            .triangle = 0,
            .kind = SparseOrdinary,
            ._reserved = 0,
            ._allocatedBySparse = 0,
        },
        .data = y,
    });
}
typedef struct {
	t_pxobject const super;
    t_chebyshev const cheby;
    t_systhread_mutex const mutex;
    double const phase;
    double const freqs;
    double const*const cache;
    double const*const coefs;
} t_sinics;

C74_HIDDEN t_class const * class = NULL;

C74_HIDDEN t_sinics const*const sinics_new(t_atom_float const value) {
    register t_sinics*this = NULL;
    if ((this = (t_sinics*const)object_alloc((t_class*const)class))) {
        z_dsp_setup((t_pxobject*const)this, 2);
        *(short*const)&this->super.z_misc |= Z_NO_INPLACE|Z_MC_INLETS;
        outlet_new(this, "signal");
        switch (systhread_mutex_new((t_systhread_mutex*const)&this->mutex, SYSTHREAD_MUTEX_NORMAL)) {
            case MAX_ERR_NONE:
                break;
            default:
                object_error((t_object*const)this, "mutex error");
                break;
        }
        *(double const**const)&this->cheby.data = NULL;
        *(double const**const)&this->coefs = (double const*const)sysmem_newptr(sizeof(double));
        *(double*const)&this->freqs = value;
        *(double*const)&this->phase = 0;
        *(double*const)this->coefs = 1;
    }
    return this;
}

C74_HIDDEN void sinics_free(t_sinics const*const this) {
    chebyshev_free((t_chebyshev*const)&this->cheby);
    if (this->coefs)
        sysmem_freeptr((void*const)this->coefs);
    systhread_mutex_free(this->mutex);
    z_dsp_free((t_pxobject*const)this);
}

C74_HIDDEN void sinics_sh(double*const y,
                          double const*const x,
                          long const length,
                          double const*const w,
                          long const count,
                          double*const z,
                          t_systhread_mutex const mutex) {
    vDSP_vsmulD(x, 1, (double const[]){2.0}, z, 1, length);
    vvcospi(y, z, (int const[]){(int const)length});
    switch (systhread_mutex_trylock(mutex)) {
        case MAX_ERR_NONE:
            vDSP_vpolyD(w+count-1, -1, y, 1, y, 1, length, count-1);
            systhread_mutex_unlock(mutex);
        default:
            break;
    }
    vvsinpi(z, z, (int const[]){(int const)length});
    vDSP_vmulD(y, 1, z, 1, y, 1, length);
}

C74_HIDDEN void sinics_dh(double*const y,
                          double const*const x,
                          double*const w,
                          long const length,
                          long const count,
                          double*const z) {
    vDSP_vsmulD(x, 1, (double const[]){2.0}, z, 1, length);
    vvcospi(y, z, (int const[]){(int const)length});
    for ( register long k = count - 1 ; 0 < k ; -- k )
        vDSP_vmaD(w+length*k, 1, y, 1, w+length*(k-1), 1, w+length*(k-1), 1, length);
    vvsinpi(z, z, (int const[]){(int const)length});
    vDSP_vmulD(w, 1, z, 1, y, 1, length);
}

C74_HIDDEN void sinics_dd(t_sinics const*const this, t_object const*const dsp64, double const*const*const i, long const ic, double*const*const o, long const oc, long const length, long const flags, void*const parameter) {
    long const count = ic - 1;
    double*const z = (double*const)(this->cache);
    double*const w = (double*const)(this->cache + length * count);
    for ( register long k = 0, K = count ; k < K ; ++ k )
        sysmem_copyptr(i[k+1], (double*const)(z+k*length), length * sizeof(double));
    chebyshev_mm(&this->cheby, w, length, z, length, length);
    double const dt = simd_precise_recip((double const)(uintptr_t const)parameter);
    vDSP_vtrapzD(*i, 1, &dt, *o, 1, length);
    vDSP_vsaddD(*o, 1, (double const[]){fma(0.5 * dt, **i, this->phase)}, *o, 1, length);
    *(double*const)&this->phase = fmod(fma(0.5 * dt, (*i)[length-1], (*o)[length-1]), 1);
    sinics_dh(*o, *o, w, length, count, z);
}

C74_HIDDEN void sinics_sd(t_sinics const*const this, t_object const*const dsp64, double const*const*const i, long const ic, double*const*const o, long const oc, long const length, long const flags, void*const parameter) {
    long const count = ic - 1;
    double*const z = (double*const)(this->cache);
    double*const w = (double*const)(this->cache + length * count);
    for ( register long k = 0, K = count ; k < K ; ++ k )
        sysmem_copyptr(i[k+1], (double*const)(z+k*length), length * sizeof(double));
    chebyshev_mm(&this->cheby, w, length, z, length, length);
    double const df = this->freqs / (uintptr_t const)parameter;
    vDSP_vrampD(&this->phase, &df, *o, 1, length);
    *(double*const)&this->phase = fmod(fma(df, length, this->phase), 1);
    sinics_dh(*o, *o, w, length, count, z);
}

C74_HIDDEN void sinics_ds(t_sinics const*const this, t_object const*const dsp64, double const*const*const i, long const ic, double*const*const o, long const oc, long const length, long const flags, void*const parameter) {
    double const dt = simd_precise_recip((double const)(uintptr_t const)parameter);
    vDSP_vtrapzD(*i, 1, &dt, *o, 1, length);
    vDSP_vsaddD(*o, 1, (double const[]){fma(0.5 * dt, **i, this->phase)}, *o, 1, length);
    *(double*const)&this->phase = fmod(fma(0.5 * dt, (*i)[length-1], (*o)[length-1]), 1);
    sinics_sh(*o, *o, length, this->coefs, sysmem_ptrsize((double*const)this->coefs) / sizeof(double), (double*const)this->cache, this->mutex);
}

C74_HIDDEN void sinics_ss(t_sinics const*const this, t_object const*const dsp64, double const*const*const i, long const ic, double*const*const o, long const oc, long const length, long const flags, void*const parameter) {
    double const df = this->freqs / (uintptr_t const)parameter;
    vDSP_vrampD(&this->phase, &df, *o, 1, length);
    *(double*const)&this->phase = fmod(fma(df, length, this->phase), 1);
    sinics_sh(*o, *o, length, this->coefs, sysmem_ptrsize((double*const)this->coefs) / sizeof(double), (double*const)this->cache, this->mutex);
}

C74_HIDDEN void sinics_dsp64(t_sinics const*const this, t_object const*const dsp64, short const*const count, double const samplerate, long const maxvectorsize, long const flags) {
    *(double*const)&this->phase = 0;
    if (count[0]&&count[1]) {
        *(double const**const)&this->cache = (double const*)sysmem_resizeptr((double*const)this->cache, 2 * this->cheby.structure.columnCount * maxvectorsize * sizeof(double));
        dsp_add64((t_object*const)dsp64, (t_object*const)this, (t_perfroutine64 const)sinics_dd, 0, (uintptr_t const)samplerate);
    } else if (count[1]) {
        *(double const**const)&this->cache = (double const*)sysmem_resizeptr((double*const)this->cache, 2 * this->cheby.structure.columnCount * maxvectorsize * sizeof(double));
        dsp_add64((t_object*const)dsp64, (t_object*const)this, (t_perfroutine64 const)sinics_sd, 0, (uintptr_t const)samplerate);
    } else if (count[0]) {
        *(double const**const)&this->cache = (double const*)sysmem_resizeptr((double*const)this->cache, maxvectorsize * sizeof(double));
        dsp_add64((t_object*const)dsp64, (t_object*const)this, (t_perfroutine64 const)sinics_ds, 0, (uintptr_t const)samplerate);
    } else {
        *(double const**const)&this->cache = (double const*)sysmem_resizeptr((double*const)this->cache, maxvectorsize * sizeof(double));
        dsp_add64((t_object*const)dsp64, (t_object*const)this, (t_perfroutine64 const)sinics_ss, 0, (uintptr_t const)samplerate);
    }
}

C74_HIDDEN long sinics_input(t_sinics const*const this, long const index, long const chans) {
    switch (index) {
        case 0:
            return chans == 1;
        case 1:
            chebyshev_init((t_chebyshev*const)&this->cheby, chans);
            return true;
    }
}

C74_HIDDEN long sinics_output(t_sinics const*const this, long const index) {
    return index == 0;
}

C74_HIDDEN void sinics_freqz(t_sinics const*const this, t_atom_float const value) {
    switch (proxy_getinlet((t_object*const)this)) {
        case 0: {
            *(double*const)&this->freqs = value;
            break;
        }
        default:
            object_error((t_object*const)this, "invalid inlet");
            break;
    }
}

C74_HIDDEN void sinics_coefs(t_sinics const*const this, t_symbol const*const symbol, short const argc, t_atom const*const argv) {
    switch (proxy_getinlet((t_object*const)this)) {
        case 1: {
            t_chebyshev chebyshev = {0};
            double*const memory = (double*const)alloca(2 * argc * sizeof(double));
            switch (atom_getdouble_array(argc, (t_atom*const)argv, argc, memory)) {
                case MAX_ERR_NONE:
                    chebyshev_init(&chebyshev, argc);
                    chebyshev_mv(&chebyshev, memory + argc, memory);
                    chebyshev_free(&chebyshev);
                    switch (systhread_mutex_lock(this->mutex)) {
                        case MAX_ERR_NONE:
                            sysmem_copyptr(memory + argc, (*(double**const)&this->coefs = (double*const)sysmem_resizeptr((double*const)this->coefs, argc * sizeof(double))), argc * sizeof(double));
                            systhread_mutex_unlock(this->mutex);
                            break;
                        default:
                            object_error((t_object*const)this, "mutex error");
                            break;
                    }
                    break;
                default:
                    object_error((t_object*const)this, "message should be array of number");
                    break;
            }
            break;
        }
        default:
            object_error((t_object*const)this, "invalid inlet");
            break;
    }
}

C74_EXPORT void ext_main(void*const r) {
    if ((class = (t_class*const)class_new("sinics~", (method const)sinics_new, (method const)sinics_free, (long const)sizeof(t_sinics), 0L, A_DEFFLOAT, 0))) {
        class_addmethod((t_class*const)class, (method const)sinics_dsp64, "dsp64", A_CANT, 0);
        class_addmethod((t_class*const)class, (method const)sinics_freqz, "float", A_FLOAT, 0);
        class_addmethod((t_class*const)class, (method const)sinics_coefs, "list", A_GIMME, 0);
        class_addmethod((t_class*const)class, (method const)sinics_input, "inputchanged", A_CANT, 0);
        class_addmethod((t_class*const)class, (method const)sinics_output, "multichanneloutputs", A_CANT, 0);
        class_dspinit((t_class*const)class);
        class_register(CLASS_BOX, (t_class*const)class);
    }
}
