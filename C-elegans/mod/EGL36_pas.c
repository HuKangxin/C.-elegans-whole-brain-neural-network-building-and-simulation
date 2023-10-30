/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__EGL36_pas
#define _nrn_initial _nrn_initial__EGL36_pas
#define nrn_cur _nrn_cur__EGL36_pas
#define _nrn_current _nrn_current__EGL36_pas
#define nrn_jacob _nrn_jacob__EGL36_pas
#define nrn_state _nrn_state__EGL36_pas
#define _net_receive _net_receive__EGL36_pas 
#define rates rates__EGL36_pas 
#define states states__EGL36_pas 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gegl36bar _p[0]
#define gegl36 _p[1]
#define m_f _p[2]
#define m_m _p[3]
#define m_s _p[4]
#define Dm_f _p[5]
#define Dm_m _p[6]
#define Dm_s _p[7]
#define ek _p[8]
#define ik _p[9]
#define v _p[10]
#define _g _p[11]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_EGL36_pas", _hoc_setdata,
 "rates_EGL36_pas", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
 static int _thread1data_inuse = 0;
static double _thread1data[4];
#define _gth 0
#define k_a k_a_EGL36_pas
 double k_a = 28.5;
#define m_stau_a m_stau_a_EGL36_pas
 double m_stau_a = 355;
#define m_mtau_a m_mtau_a_EGL36_pas
 double m_mtau_a = 63;
#define m_ftau_a m_ftau_a_EGL36_pas
 double m_ftau_a = 13;
#define m_stau_EGL36_pas _thread1data[0]
#define m_stau _thread[_gth]._pval[0]
#define m_mtau_EGL36_pas _thread1data[1]
#define m_mtau _thread[_gth]._pval[1]
#define m_ftau_EGL36_pas _thread1data[2]
#define m_ftau _thread[_gth]._pval[2]
#define m_V m_V_EGL36_pas
 double m_V = 63;
#define minf_EGL36_pas _thread1data[3]
#define minf _thread[_gth]._pval[3]
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "m_V_EGL36_pas", "mV",
 "k_a_EGL36_pas", "mV",
 "m_stau_a_EGL36_pas", "ms",
 "m_mtau_a_EGL36_pas", "ms",
 "m_ftau_a_EGL36_pas", "ms",
 "m_ftau_EGL36_pas", "ms",
 "m_mtau_EGL36_pas", "ms",
 "m_stau_EGL36_pas", "ms",
 "gegl36bar_EGL36_pas", "nS/cm2",
 "gegl36_EGL36_pas", "nS/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double m_s0 = 0;
 static double m_m0 = 0;
 static double m_f0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "m_V_EGL36_pas", &m_V_EGL36_pas,
 "k_a_EGL36_pas", &k_a_EGL36_pas,
 "m_stau_a_EGL36_pas", &m_stau_a_EGL36_pas,
 "m_mtau_a_EGL36_pas", &m_mtau_a_EGL36_pas,
 "m_ftau_a_EGL36_pas", &m_ftau_a_EGL36_pas,
 "minf_EGL36_pas", &minf_EGL36_pas,
 "m_ftau_EGL36_pas", &m_ftau_EGL36_pas,
 "m_mtau_EGL36_pas", &m_mtau_EGL36_pas,
 "m_stau_EGL36_pas", &m_stau_EGL36_pas,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"EGL36_pas",
 "gegl36bar_EGL36_pas",
 0,
 "gegl36_EGL36_pas",
 0,
 "m_f_EGL36_pas",
 "m_m_EGL36_pas",
 "m_s_EGL36_pas",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 12, _prop);
 	/*initialize range parameters*/
 	gegl36bar = 1.3;
 	_prop->param = _p;
 	_prop->param_size = 12;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _EGL36_pas_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 2);
  _extcall_thread = (Datum*)ecalloc(1, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
  _thread1data_inuse = 0;
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 12, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 EGL36_pas D:/download/yanjjiusheng/hukangxin/hukangxin/mod/EGL36_pas.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "EGL36_pas.mod Potassium ion channel of HH model in C.elegans";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[3], _dlist1[3];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Dm_f = ( minf - m_f ) / m_ftau ;
   Dm_m = ( minf - m_m ) / m_mtau ;
   Dm_s = ( minf - m_s ) / m_stau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargscomma_ v ) ;
 Dm_f = Dm_f  / (1. - dt*( ( ( ( - 1.0 ) ) ) / m_ftau )) ;
 Dm_m = Dm_m  / (1. - dt*( ( ( ( - 1.0 ) ) ) / m_mtau )) ;
 Dm_s = Dm_s  / (1. - dt*( ( ( ( - 1.0 ) ) ) / m_stau )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargscomma_ v ) ;
    m_f = m_f + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / m_ftau)))*(- ( ( ( minf ) ) / m_ftau ) / ( ( ( ( - 1.0 ) ) ) / m_ftau ) - m_f) ;
    m_m = m_m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / m_mtau)))*(- ( ( ( minf ) ) / m_mtau ) / ( ( ( ( - 1.0 ) ) ) / m_mtau ) - m_m) ;
    m_s = m_s + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / m_stau)))*(- ( ( ( minf ) ) / m_stau ) / ( ( ( ( - 1.0 ) ) ) / m_stau ) - m_s) ;
   }
  return 0;
}
 
static int  rates ( _threadargsprotocomma_ double _lv ) {
   double _lq10 ;
  _lq10 = pow( 3.0 , ( ( celsius - 6.3 ) / 10.0 ) ) ;
   minf = 1.0 / ( 1.0 + exp ( ( m_V - _lv ) / k_a ) ) ;
   m_ftau = m_ftau_a ;
   m_mtau = m_mtau_a ;
   m_stau = m_stau_a ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_mem_init(Datum* _thread) {
  if (_thread1data_inuse) {_thread[_gth]._pval = (double*)ecalloc(4, sizeof(double));
 }else{
 _thread[_gth]._pval = _thread1data; _thread1data_inuse = 1;
 }
 }
 
static void _thread_cleanup(Datum* _thread) {
  if (_thread[_gth]._pval == _thread1data) {
   _thread1data_inuse = 0;
  }else{
   free((void*)_thread[_gth]._pval);
  }
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  m_s = m_s0;
  m_m = m_m0;
  m_f = m_f0;
 {
   rates ( _threadargscomma_ v ) ;
   m_f = 0.0 ;
   m_m = 0.0 ;
   m_s = 0.0 ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ek = _ion_ek;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gegl36 = gegl36bar * ( 0.33 * m_f + 0.36 * m_m + 0.39 * m_s ) ;
   ik = gegl36 * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ek = _ion_ek;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ek = _ion_ek;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m_f) - _p;  _dlist1[0] = &(Dm_f) - _p;
 _slist1[1] = &(m_m) - _p;  _dlist1[1] = &(Dm_m) - _p;
 _slist1[2] = &(m_s) - _p;  _dlist1[2] = &(Dm_s) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "EGL36_pas.mod";
static const char* nmodl_file_text = 
  "TITLE EGL36_pas.mod Potassium ion channel of HH model in C.elegans\n"
  "\n"
  "COMMENT\n"
  "1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.\n"
  "2. Refer to \"Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD\", Plos One, 14(7), Jul. 2019.\n"
  "3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.\n"
  "4. gegl36bar, ek, celsius\n"
  "ENDCOMMENT\n"
  "\n"
  "? interface\n"
  "NEURON {\n"
  "        SUFFIX EGL36_pas\n"
  "        USEION k READ ek WRITE ik    : ek is reversal voltage\n"
  "        RANGE  gegl36bar, gegl36     : gegl36bar is max gegl36\n"
  "        GLOBAL minf, m_V, k_a, m_ftau, m_mtau, m_stau, m_ftau_a, m_mtau_a, m_stau_a\n"
  "	THREADSAFE : assigned GLOBALs will be per thread\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "        (pA) = (picoamp)\n"
  "        (mV) = (millivolt)\n"
  "	(nS) = (nanosiemens)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "        gegl36bar  = 1.3 (nS/cm2)   : 1.3 for RMD\n"
  "    : minf\n"
  "	m_V = 63.0 (mV)\n"
  "        k_a = 28.5 (mV)\n"
  "    : mtau\n"
  "        m_stau_a    =  355.0 (ms)\n"
  "        m_mtau_a    =  63.0 (ms)\n"
  "        m_ftau_a    =  13.0 (ms)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "        m_f m_m  m_s\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "        v (mV)\n"
  "        celsius (degC)\n"
  "        ek (mV)\n"
  "        gegl36 (nS/cm2)\n"
  "	ik (pA/cm2)\n"
  "        minf\n"
  "	m_ftau (ms)\n"
  "	m_mtau (ms)\n"
  "	m_stau (ms)\n"
  "}\n"
  "\n"
  "? currents\n"
  "BREAKPOINT {\n"
  "        SOLVE states METHOD cnexp\n"
  "        gegl36 = gegl36bar*(0.33*m_f+0.36*m_m+0.39*m_s)\n"
  "	ik = gegl36*(v -ek)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	rates(v)\n"
  "	m_f = 0\n"
  "	m_m = 0\n"
  "	m_s = 0\n"
  "}\n"
  "\n"
  "? states\n"
  "DERIVATIVE states {\n"
  "        rates(v)\n"
  "        m_f' = (minf-m_f)/m_ftau\n"
  "        m_m' = (minf-m_m)/m_mtau\n"
  "        m_s' = (minf-m_s)/m_stau\n"
  "}\n"
  "\n"
  ":LOCAL q10\n"
  "\n"
  "? rates\n"
  "PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.\n"
  "                          :Call once from HOC to initialize inf at resting v.\n"
  "        LOCAL  q10\n"
  "        : `TABLE minf` will be replaced with `:TABLE minf` if CoreNEURON enabled)\n"
  "        : TABLE minf, m_ftau, m_mtau, m_stau DEPEND celsius FROM -100 TO 100 WITH 200\n"
  "\n"
  "UNITSOFF\n"
  "        q10 = 3^((celsius - 6.3)/10)\n"
  "\n"
  "                :\"m\" EGL-36 type potassium activation system\n"
  "		minf = 1/(1+exp((m_V-v)/k_a))\n"
  "        m_ftau = m_ftau_a   : multiplied by q10, if consider temperature\n"
  "		m_mtau = m_mtau_a   : multiplied by q10, if consider temperature\n"
  "		m_stau = m_stau_a   : multiplied by q10, if consider temperature\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif
