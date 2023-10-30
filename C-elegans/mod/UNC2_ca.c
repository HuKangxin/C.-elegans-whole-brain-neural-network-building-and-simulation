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
 
#define nrn_init _nrn_init__UNC2_ca
#define _nrn_initial _nrn_initial__UNC2_ca
#define nrn_cur _nrn_cur__UNC2_ca
#define _nrn_current _nrn_current__UNC2_ca
#define nrn_jacob _nrn_jacob__UNC2_ca
#define nrn_state _nrn_state__UNC2_ca
#define _net_receive _net_receive__UNC2_ca 
#define rates rates__UNC2_ca 
#define states states__UNC2_ca 
 
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
#define gunc2bar _p[0]
#define gunc2 _p[1]
#define minf _p[2]
#define hinf _p[3]
#define mtau _p[4]
#define htau _p[5]
#define m _p[6]
#define h _p[7]
#define Dm _p[8]
#define Dh _p[9]
#define eca _p[10]
#define ica _p[11]
#define v _p[12]
#define _g _p[13]
#define _ion_eca	*_ppvar[0]._pval
#define _ion_ica	*_ppvar[1]._pval
#define _ion_dicadv	*_ppvar[2]._pval
 
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
 "setdata_UNC2_ca", _hoc_setdata,
 "rates_UNC2_ca", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
#define h_f h_f_UNC2_ca
 double h_f = -3.6;
#define h_e h_e_UNC2_ca
 double h_e = -6.1;
#define h_d h_d_UNC2_ca
 double h_d = 122.6;
#define h_c h_c_UNC2_ca
 double h_c = -3.5;
#define h_b h_b_UNC2_ca
 double h_b = 22.9;
#define h_a h_a_UNC2_ca
 double h_a = 142.5;
#define h_V h_V_UNC2_ca
 double h_V = -77.5;
#define k_i k_i_UNC2_ca
 double k_i = 5.6;
#define k_a k_a_UNC2_ca
 double k_a = 4;
#define m_e m_e_UNC2_ca
 double m_e = 0.1;
#define m_d m_d_UNC2_ca
 double m_d = 15.4;
#define m_c m_c_UNC2_ca
 double m_c = 9.1;
#define m_b m_b_UNC2_ca
 double m_b = -38.2;
#define m_a m_a_UNC2_ca
 double m_a = 1.5;
#define m_V m_V_UNC2_ca
 double m_V = -37.2;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "m_V_UNC2_ca", "mV",
 "k_a_UNC2_ca", "mV",
 "h_V_UNC2_ca", "mV",
 "k_i_UNC2_ca", "mV",
 "m_a_UNC2_ca", "ms",
 "m_b_UNC2_ca", "mV",
 "m_c_UNC2_ca", "mV",
 "m_d_UNC2_ca", "mV",
 "m_e_UNC2_ca", "ms",
 "h_a_UNC2_ca", "ms",
 "h_b_UNC2_ca", "mV",
 "h_c_UNC2_ca", "mV",
 "h_d_UNC2_ca", "ms",
 "h_e_UNC2_ca", "mV",
 "h_f_UNC2_ca", "mV",
 "gunc2bar_UNC2_ca", "nS/cm2",
 "gunc2_UNC2_ca", "nS/cm2",
 "mtau_UNC2_ca", "ms",
 "htau_UNC2_ca", "ms",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "m_V_UNC2_ca", &m_V_UNC2_ca,
 "k_a_UNC2_ca", &k_a_UNC2_ca,
 "h_V_UNC2_ca", &h_V_UNC2_ca,
 "k_i_UNC2_ca", &k_i_UNC2_ca,
 "m_a_UNC2_ca", &m_a_UNC2_ca,
 "m_b_UNC2_ca", &m_b_UNC2_ca,
 "m_c_UNC2_ca", &m_c_UNC2_ca,
 "m_d_UNC2_ca", &m_d_UNC2_ca,
 "m_e_UNC2_ca", &m_e_UNC2_ca,
 "h_a_UNC2_ca", &h_a_UNC2_ca,
 "h_b_UNC2_ca", &h_b_UNC2_ca,
 "h_c_UNC2_ca", &h_c_UNC2_ca,
 "h_d_UNC2_ca", &h_d_UNC2_ca,
 "h_e_UNC2_ca", &h_e_UNC2_ca,
 "h_f_UNC2_ca", &h_f_UNC2_ca,
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
"UNC2_ca",
 "gunc2bar_UNC2_ca",
 0,
 "gunc2_UNC2_ca",
 "minf_UNC2_ca",
 "hinf_UNC2_ca",
 "mtau_UNC2_ca",
 "htau_UNC2_ca",
 0,
 "m_UNC2_ca",
 "h_UNC2_ca",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 14, _prop);
 	/*initialize range parameters*/
 	gunc2bar = 1;
 	_prop->param = _p;
 	_prop->param_size = 14;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* eca */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _UNC2_ca_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 14, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 UNC2_ca D:/download/yanjjiusheng/hukangxin/hukangxin/mod/UNC2_ca.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "UNC2_ca.mod Potassium ion channel of HH model in C.elegans";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0 ) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0 ) ) ) / htau ) - h) ;
   }
  return 0;
}
 
static int  rates ( _threadargsprotocomma_ double _lv ) {
   double _lq10 ;
  _lq10 = pow( 3.0 , ( ( celsius - 6.3 ) / 10.0 ) ) ;
   minf = 1.0 / ( 1.0 + exp ( ( m_V - _lv ) / k_a ) ) ;
   mtau = m_a / ( exp ( ( m_b - _lv ) / m_c ) + exp ( ( _lv - m_b ) / m_d ) ) + m_e ;
   hinf = 1.0 / ( 1.0 + exp ( ( _lv - h_V ) / k_i ) ) ;
   htau = h_a / ( 1.0 + exp ( ( h_b - _lv ) / h_c ) ) + h_d / ( 1.0 + exp ( ( _lv - h_e ) / h_f ) ) ;
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
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  eca = _ion_eca;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
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
  eca = _ion_eca;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
 {
   rates ( _threadargscomma_ v ) ;
   m = 0.0 ;
   h = 1.0 ;
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
  eca = _ion_eca;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gunc2 = gunc2bar * m * h ;
   ica = gunc2 * ( v - 45.0 ) ;
   }
 _current += ica;

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
  eca = _ion_eca;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
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
  eca = _ion_eca;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "UNC2_ca.mod";
static const char* nmodl_file_text = 
  "TITLE UNC2_ca.mod Potassium ion channel of HH model in C.elegans\n"
  "\n"
  "COMMENT\n"
  "1. Created by Yu Zhang at Tsinghua University, on Sep. 2021.\n"
  "2. Refer to \"Biophysical modeling of C. elegans neurons: Single ion currents and whole-cell dynamics of AWC on and RMD\", Plos One, 14(7), Jul. 2019.\n"
  "3. Membrane voltage is in absolute mV and has been reversed in polarity from the original HH convention.\n"
  "4. gunc2bar, eca, celsius\n"
  "ENDCOMMENT\n"
  "\n"
  "? interface\n"
  "NEURON {\n"
  "        SUFFIX UNC2_ca\n"
  "        USEION ca READ eca WRITE ica    	  : eca is reversal voltage\n"
  "        RANGE  gunc2bar, gunc2, minf, mtau, hinf, htau\n"
  "		GLOBAL m_V, k_a, h_V, k_i, m_a, m_b, m_c, m_d, m_e, h_a, h_b, h_c, h_d, h_e, h_f\n"
  "		THREADSAFE : assigned GLOBALs will be per thread\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "        (pA) = (picoamp)\n"
  "        (mV) = (millivolt)\n"
  "	    (nS) = (nanosiemens)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "        gunc2bar  = 1.0 (nS/cm2)   : 1.0 for AWC_ON, 0.9 for RMD\n"
  "    : minf\n"
  "		m_V    = -37.2 (mV)        : -37.2 as an exchange\n"
  "        k_a    = 4.0 (mV)\n"
  "    : hinf\n"
  "        h_V    = -77.5 (mV)        : -77.5 as an exchange\n"
  "        k_i    = 5.6 (mV)\n"
  "    : mtau\n"
  "        m_a    = 1.5 (ms)\n"
  "        m_b    = -38.2 (mV)       : 3.0 as an exchange\n"
  "        m_c    = 9.1 (mV)\n"
  "        m_d    = 15.4 (mV)         : 0.1 as an exchange\n"
  "        m_e    = 0.1 (ms)         : 0.1 as an exchange\n"
  "    : htau\n"
  "        h_a    = 142.5 (ms)        : 142.5 as an exchange\n"
  "        h_b    = 22.9 (mV)        : 22.9 as an exchange\n"
  "        h_c    = -3.5 (mV)\n"
  "        h_d    = 122.6 (ms)        : 122.6 as an exchange\n"
  "        h_e    = -6.1 (mV)        : -6.1 as an exchange\n"
  "        h_f    = -3.6  (mV)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "        m h\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "        v (mV)\n"
  "        celsius (degC)\n"
  "        eca (mV)\n"
  "        gunc2 (nS/cm2)\n"
  "		ica (pA/cm2)\n"
  "        minf hinf\n"
  "	    mtau (ms) htau (ms)\n"
  "}\n"
  "\n"
  "? currents\n"
  "BREAKPOINT {\n"
  "        SOLVE states METHOD cnexp\n"
  "        gunc2 = gunc2bar*m*h\n"
  "	    ica = gunc2*(v - 45)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	rates(v)\n"
  "	m = 0  : minf\n"
  "	h = 1  : hinf\n"
  "}\n"
  "\n"
  "? states\n"
  "DERIVATIVE states {\n"
  "        rates(v)\n"
  "        m' = (minf-m)/mtau\n"
  "        h' = (hinf-h)/htau\n"
  "}\n"
  "\n"
  ":LOCAL q10\n"
  "\n"
  "? rates\n"
  "PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.\n"
  "                          :Call once from HOC to initialize inf at resting v.\n"
  "        LOCAL  q10\n"
  "        : `TABLE minf` will be replaced with `:TABLE minf` if CoreNEURON enabled)\n"
  "        : TABLE minf, mtau, hinf, htau DEPEND celsius FROM -100 TO 100 WITH 200\n"
  "\n"
  "UNITSOFF\n"
  "        q10 = 3^((celsius - 6.3)/10)\n"
  "\n"
  "                :\"m\" UNC-2 type potassium activation system\n"
  "		minf = 1/(1+exp((m_V-v)/k_a))\n"
  "        mtau = m_a/(exp((m_b-v)/m_c)+exp((v-m_b)/m_d))+m_e   : multiplied by q10, if consider temperature\n"
  "\n"
  "                :\"h\" UNC-2 type potassium inactivation system\n"
  "		hinf = 1/(1+exp((v-h_V)/k_i))\n"
  "		htau = h_a/(1+exp((h_b-v)/h_c))+h_d/(1+exp((v-h_e)/h_f))   : multiplied by q10, if consider temperature\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif
