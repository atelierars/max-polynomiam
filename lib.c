#include"ext.h"
#include"ext_obex.h"
#include"z_dsp.h"
#include<simd/simd.h>
#include<Accelerate/Accelerate.h>
typedef SparseMatrix_Double t_sparsematrix;
C74_HIDDEN void sparse_chebyshev1(t_sparsematrix*const matrix, long const count) {
	register long const nz = ( ( count / 2 ) + 1 ) * ( ( count + 1 ) / 2 );
	void * const memory = sysmem_resizeptrclear(matrix->data, ( count + 1 ) * sizeof(long const) + nz * (sizeof(double const) + sizeof(int const)));
	matrix->structure.rowCount = (int const)count;
	matrix->structure.columnCount = (int const)count;
	matrix->structure.rowIndices = (int*const)(memory + nz * (sizeof(double const)));
	matrix->structure.columnStarts = (long*const)(memory + nz * (sizeof(double const) + sizeof(int const)));
	matrix->structure.attributes = (SparseAttributes_t const) {
		.transpose = 0,
		.triangle = SparseUpperTriangle,
		.kind = SparseTriangular,
		._reserved = 0,
		._allocatedBySparse = 0,
	};
	matrix->structure.blockSize = 1;
	matrix->data = (double*const)memory;
	switch (count) {
		case 0:
		case 1:
			matrix->data[0] = 1;
			matrix->structure.rowIndices[0] = 0;
			matrix->structure.columnStarts[0] = 0;
			matrix->structure.columnStarts[1] = 1;
			break;
		default:
			vDSP_vclrD(matrix->data, 1, nz);
			matrix->data[0] = 1;
			matrix->data[1] = 1;
			matrix->structure.rowIndices[0] = 0;
			matrix->structure.rowIndices[1] = 1;
			matrix->structure.columnStarts[0] = 0;
			matrix->structure.columnStarts[1] = 1;
			matrix->structure.columnStarts[2] = 2;
			for ( register unsigned long c = 2, C = count ; c < C ; ++ c ) {
				register unsigned long const j = matrix->structure.columnStarts[c-2];
				register unsigned long const k = matrix->structure.columnStarts[c-1];
				register unsigned long const l = matrix->structure.columnStarts[c-0];
				register unsigned long const p = c % 2;
				register unsigned long const q = c / 2 + 1;
				for ( register unsigned long r = 0, R = q ; r < R ; ++ r )
					matrix->structure.rowIndices[l+r] = (int const)(2 * r + p);
				vDSP_vsmulD(matrix->data + k, 1, (double const[]){2.0}, matrix->data + l + 1 - p, 1, l - k);
				vDSP_vsubD(matrix->data + j, 1, matrix->data + l, 1, matrix->data + l, 1, k - j);
				matrix->structure.columnStarts[c+1] = l + q;
			}
			break;
	}
}
C74_HIDDEN void sparse_destroy(t_sparsematrix*const matrix) {
	sysmem_freeptr(matrix->data);
	matrix->data = NULL;
}
// Col-Major Vector <- SparseMatrix * Col-Major Vector
C74_HIDDEN void sparse_gemv(t_sparsematrix const*const matrix, double*const y, double*const x) {
	SparseMultiply(*matrix,
				   (DenseVector_Double const) { .data = x, .count = matrix->structure.columnCount, },
				   (DenseVector_Double const) { .data = y, .count = matrix->structure.rowCount, });
}
// Col-Major Matrix <- SparseMatrix * Col-Major Matrix
C74_HIDDEN void sparse_gemm_cc(t_sparsematrix const*const matrix, double*const y, long const ldy, double*const x, long const ldx, long const length) {
	SparseMultiply(*matrix,
				   (DenseMatrix_Double const) {
		.rowCount = matrix->structure.columnCount,
		.columnCount = (int const)length,
		.columnStride = (int const)ldx,
		.attributes = (SparseAttributes_t const) {
			.transpose = 0,
			.triangle = 0,
			.kind = SparseOrdinary,
			._reserved = 0,
			._allocatedBySparse = 0,
		},
		.data = x,
	},
				   (DenseMatrix_Double const) {
		.rowCount = matrix->structure.rowCount,
		.columnCount = (int const)length,
		.columnStride = (int const)ldy,
		.attributes = (SparseAttributes_t const) {
			.transpose = 0,
			.triangle = 0,
			.kind = SparseOrdinary,
			._reserved = 0,
			._allocatedBySparse = 0,
		},
		.data = y,
	});
}
// Row-Major Matrix <- SparseMatrix * Row-Major Matrix
C74_HIDDEN void sparse_gemm_rr(t_sparsematrix const*const matrix, double*const y, long const ldy, double*const x, long const ldx, long const length) {
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
	t_sparsematrix const cheby;
	t_systhread_mutex const mutex;
	uintptr_t const split;
	double*const cache;
	double*const coefs;
} t_polynomiam;
C74_HIDDEN t_class const * class = NULL;
C74_HIDDEN t_polynomiam const*const new(void) {
	register t_polynomiam * const object = (t_polynomiam*const)object_alloc((t_class*const)class);
	if ( object ) {
		z_dsp_setup((t_pxobject*const)object, 2);
		*(short*const)&object->super.z_misc |= Z_MC_INLETS;//|Z_NO_INPLACE;
		outlet_new(object, "signal");
		switch (systhread_mutex_new((t_systhread_mutex*const)&object->mutex, SYSTHREAD_MUTEX_NORMAL)) {
			case MAX_ERR_NONE:
				break;
			default:
				object_error((t_object*const)object, "mutex error");
				break;
		}
		*(double const**const)&object->cheby.data = NULL;
		*(double const**const)&object->coefs = (double const*const)sysmem_newptr(sizeof(double const));
		*(double const**const)&object->cache = (double const*const)sysmem_newptr(0);
		*(double*const)object->coefs = 1;
	}
	return object;
}
C74_HIDDEN void del(t_polynomiam const*const this) {
	z_dsp_free((t_pxobject*const)this);
	systhread_mutex_free(this->mutex);
	sysmem_freeptr((void*const)this->coefs);
	sysmem_freeptr((void*const)this->cache);
	sparse_destroy((t_sparsematrix*const)&this->cheby);
}
C74_HIDDEN void clr(t_polynomiam const*const this, t_object const*const dsp64, double const*const*const ins, long const numins, double*const*const outs, long const numouts, long const length, long const flags, void*const parameter) {
	for ( register long k = 0, K = numouts ; k < K ; ++ k )
		vDSP_vclrD(outs[k], 1, length);
}
C74_HIDDEN void fix(t_polynomiam const*const this, t_object const*const dsp64, double const*const*const ins, long const numins, double*const*const outs, long const numouts, long const length, long const flags, void*const parameter) {
	switch (systhread_mutex_trylock(this->mutex)) {
		case MAX_ERR_NONE:
			for ( register unsigned long k = 0, K = simd_reduce_min((simd_ulong2 const){numins, numouts}), count = sysmem_ptrsize((void*const)this->coefs) / sizeof(double const) ; k < K ; ++ k )
				vDSP_vpolyD(this->coefs+count-1, -1, ins[k], 1, outs[k], 1, length, count-1);
			systhread_mutex_unlock(this->mutex);
			break;
	}
}
C74_HIDDEN void dyn(t_polynomiam const*const this, t_object const*const dsp64, double const*const*const ins, long const numins, double*const*const outs, long const numouts, long const length, long const flags, void*const parameter) {
	register uintptr_t const ld = (uintptr_t const)parameter;
	register uintptr_t const np = this->cheby.structure.columnCount; // number of harmonics
	register uintptr_t const sp = this->split;
	register double * const z = this->cache + np * ld;
	register double * const w = this->cache;
	assert(numins == np + sp);
	assert(length <= ld);
	for ( register uintptr_t k = 0, K = np ; k < K ; ++ k )
		sysmem_copyptr(ins[sp+k], z + ld * k, length * sizeof(double const));
	sparse_gemm_rr(&this->cheby, w, ld, z, ld, length);
	sysmem_copyptr(*ins, z, length * sizeof(double const));
	for ( register uintptr_t k = 1, K = np ; k < K ; ++ k, vDSP_vmulD(*ins, 1, z, 1, z, 1, length) )
		vDSP_vmaD(z, 1, w + ld * k, 1, w, 1, w, 1, length);
	sysmem_copyptr(w, *outs, length * sizeof(double const));
}
C74_HIDDEN void dsp64(t_polynomiam const*const this, t_object const*const dsp64, short const*const count, double const samplerate, long const maxvectorsize, long const flags) {
	register long const length = this->cheby.structure.columnCount * 2 * maxvectorsize;
	*(double const**const)&this->cache = (double*const)sysmem_resizeptr((void*const)this->cache, length * sizeof(double const));
	dsp_add64((t_object*const)dsp64, (t_object*const)this, (t_perfroutine64 const)(count[0] ? count[1] ? dyn : fix : clr), 0, (uintptr_t const)maxvectorsize);
}
C74_HIDDEN void list(t_polynomiam const*const this, t_symbol const*const symbol, short const argc, t_atom const*const argv) {
	t_sparsematrix chebyshev = {0};
	double*const memory = (double*const)alloca(2 * argc * sizeof(double const));
	switch (proxy_getinlet((t_object*const)this)) {
		case 1:
			switch (atom_getdouble_array(argc, (t_atom*const)argv, argc, memory)) {
				case MAX_ERR_NONE:
					sparse_chebyshev1(&chebyshev, argc);
					sparse_gemv(&chebyshev, memory + argc, memory);
					sparse_destroy(&chebyshev);
					switch (systhread_mutex_lock(this->mutex)) {
						case MAX_ERR_NONE:
							sysmem_copyptr(memory + argc, *(double**const)&this->coefs = (double*const)sysmem_resizeptr((double*const)this->coefs, argc * sizeof(double const)), argc * sizeof(double const));
							systhread_mutex_unlock(this->mutex);
							break;
						default:
							object_error((t_object*const)this, "mutex error");
							break;
					}
					break;
				default:
					object_error((t_object*const)this, "coefficient should be array of number");
					break;
			}
			break;
		default:
			object_error((t_object*const)this, "invalid inlet");
			break;
	}
}
C74_HIDDEN long input(t_polynomiam const*const this, long const index, long const chans) {
	switch (index) {
		case 0:
			*(uintptr_t*const)&this->split = chans;
			return chans == 1;
		case 1:
			sparse_chebyshev1((t_sparsematrix*const)&this->cheby, chans);
			return true;
	}
}
C74_HIDDEN long output(t_polynomiam const*const this, long const index) {
	return!index;
}
C74_EXPORT void ext_main(void*const r) {
	if (!class) {
		t_class * const object = (t_class*const)class_new("polynomiam~", (method const)new, (method const)del, (long const)sizeof(t_polynomiam const), 0L, 0L);
		class_addmethod((t_class*const)object, (method const)dsp64, "dsp64", A_CANT, 0);
		class_addmethod((t_class*const)object, (method const)list, "list", A_GIMME, 0);
		class_addmethod((t_class*const)object, (method const)input, "inputchanged", A_CANT, 0);
		class_addmethod((t_class*const)object, (method const)output, "multichanneloutputs", A_CANT, 0);
		class_dspinit((t_class*const)object);
		class_register(CLASS_BOX, (t_class*const)object);
		class = object;
	}
}


