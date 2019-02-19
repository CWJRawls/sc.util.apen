/**
	@file
	sc.util.apen - an object to compute the approximate entropy of a time seires of data
	Connor Rawls - cwrawls@asu.edu

    Copyright Synthesis Center, Arizona State University, 2018
 
	@ingroup    analysis-utilities
*/

#include "ext.h"							// standard Max include, always required
#include "ext_obex.h"						// required for new style Max object

////////////////////////// object struct
typedef struct _sc_util_apen
{
	t_object	            ob;
    long                    series_length;              //the current size of the array. Must be >= 1 <= series_max_length
    long                    series_max_length;          //the maximum size of the array, will be assigned a default value
    long                    series_vector_size;         //the size of the vector held at each point in the series
    double                  similarity;                 //the thresholding factor when considering the similarity between patterns
    long                    pattern_length;             //the number of points in the series considered in a single pattern
    long                    calc_on_input;              //flag to determine if ApEn should be calculated whenever new input is received
    long                    hold_size_warning;          //flag to determine if ApEn should print to the console when there is insufficient data to compute
    double*                 test_value;                 //holds data series. Will replace with Eigen Array/Matrix when moving to N-D vectors
	void		            *out;                       //outlet
    void*                   out2;                       //dumpout
} t_sc_util_apen;

///////////////////////// function prototypes
//// standard set
void *sc_util_apen_new(t_symbol *s, long argc, t_atom *argv);
void sc_util_apen_free(t_sc_util_apen *x);
void sc_util_apen_assist(t_sc_util_apen *x, void *b, long m, long a, char *s);

//// For handling input
void sc_util_apen_bang(t_sc_util_apen *x);                                          //Regardless of calc_on_input calculates ApEn if enough data is available
void sc_util_apen_anything(t_sc_util_apen *x, t_symbol *s, long ac, t_atom *av);    //Checks data content and size, adds to series if data is in the correct format and length
void sc_util_apen_clear(t_sc_util_apen *x);                                         //empties list of data
void sc_util_apen_set_series_length(t_sc_util_apen *x, void *attr, long argc, t_atom *argv);                     //sets the maximum length of the series
void sc_util_apen_set_vector_size(t_sc_util_apen *x, void *attr, long argc, t_atom *argv);                       //sets the size of the data vector and clears the list
void sc_util_apen_pattern_length(t_sc_util_apen *x, void *attr, long argc, t_atom *argv);                        //sets the size of the pattern to be computed
void sc_util_apen_similarity(t_sc_util_apen *x, void *attr, long argc, t_atom *argv);                           //sets the threshold for pattern similarity
void sc_util_apen_calc_on_input(t_sc_util_apen *x, void *attr, long argc, t_atom *argv);                         //sets whether or not to attempt calculating ApEn when a new data point is received
void sc_util_apen_hold_size_warning(t_sc_util_apen *x, void *attr, long argc, t_atom *argv);                     //sets flag for showing insufficient data warnings

t_max_err sc_util_apen_notify(t_sc_util_apen *x, t_symbol *s, t_symbol *msg, void *sender, void *data);

//Attribute Getters
void sc_util_apen_get_similarity(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv);
void sc_util_apen_get_coi(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv);
void sc_util_apen_get_series_length(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv);
void sc_util_apen_get_cur_size(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv);
void sc_util_apen_set_cur_size(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv);       //dummy function to prevent attribute being set
void sc_util_apen_get_vector_size(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv);
void sc_util_apen_get_pattern_length(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv);
void sc_util_apen_get_size_warning(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv);


void sc_util_apen_dump(t_sc_util_apen *x); //Get a list of stored values out the right outlet

void sc_util_apen_calculate(t_sc_util_apen *x); //function to actually calculate Approximate Entropy

double sc_util_apen_maxdist(double* d0, double* d1, long l, double r); //get the maximum distance between pattern components
void sc_util_apen_getstate(t_sc_util_apen* x); //output all values through the dumpout

//Functions for inputting new data
void sc_util_apen_int(t_sc_util_apen *x, long n);
void sc_util_apen_float(t_sc_util_apen *x, double f);
void sc_util_apen_list(t_sc_util_apen *x, t_symbol* a, long argc, t_atom *argv);

//////////////////////// global class pointer variable
void *sc_util_apen_class;


void ext_main(void *r)
{
	t_class *c;

	c = class_new("sc.apen", (method)sc_util_apen_new, (method)sc_util_apen_free, (long)sizeof(t_sc_util_apen),
				  0L /* leave NULL!! */, A_GIMME, 0);

	class_addmethod(c, (method)sc_util_apen_bang,			    "bang",                             0);
    class_addmethod(c, (method)sc_util_apen_clear,              "clear",                            0);
    class_addmethod(c, (method)sc_util_apen_dump,               "dump",                             0);
    //class_addmethod(c, (method)sc_util_apen_hold_size_warning,  "size_warning",         A_LONG,     0);
    class_addmethod(c, (method)sc_util_apen_anything,           "anything",             A_GIMME,    0);
    class_addmethod(c, (method)sc_util_apen_int,                "int",                  A_LONG,     0);
    class_addmethod(c, (method)sc_util_apen_notify,             "notify",               A_CANT,     0);
    class_addmethod(c, (method)sc_util_apen_float,              "float",                A_FLOAT,    0);
    class_addmethod(c, (method)sc_util_apen_getstate,           "getstate",                         0);
    class_addmethod(c, (method)sc_util_apen_list,               "list",                 A_GIMME,    0);
    
    //Symbol versions of attributes we want to be callable from the patcher
    CLASS_ATTR_LONG(c, "series_length",          0,                      t_sc_util_apen , series_max_length);
    CLASS_ATTR_ACCESSORS(c, "series_length", sc_util_apen_get_series_length, sc_util_apen_set_series_length);
    
    CLASS_ATTR_LONG(c, "current_size",           ATTR_SET_OPAQUE,        t_sc_util_apen, series_length);
    CLASS_ATTR_ACCESSORS(c, "current_size", sc_util_apen_get_cur_size, sc_util_apen_set_cur_size);
    
    CLASS_ATTR_LONG(c, "vector_size",            0,                      t_sc_util_apen, series_vector_size);
    CLASS_ATTR_ACCESSORS(c, "vector_size", sc_util_apen_get_vector_size,sc_util_apen_set_vector_size);
    
    CLASS_ATTR_LONG(c, "pattern_length",         0,                      t_sc_util_apen, pattern_length);
    CLASS_ATTR_ACCESSORS(c, "pattern_length", sc_util_apen_get_pattern_length, sc_util_apen_pattern_length);
    
    CLASS_ATTR_DOUBLE(c, "similarity",             0,                      t_sc_util_apen, similarity);
    CLASS_ATTR_ACCESSORS(c, "similarity",        sc_util_apen_get_similarity,       sc_util_apen_similarity);
    
    CLASS_ATTR_LONG(c, "calculate_on_input",     0,                      t_sc_util_apen, calc_on_input);
    CLASS_ATTR_STYLE(c, "calculate_on_input",   0,                       "onoff");
    CLASS_ATTR_ACCESSORS(c, "calculate_on_input", sc_util_apen_get_coi, sc_util_apen_calc_on_input);
    
    CLASS_ATTR_LONG(c, "size_warning",           0,                      t_sc_util_apen, hold_size_warning);
    CLASS_ATTR_STYLE(c, "size_warning", 0, "onoff");
    CLASS_ATTR_ACCESSORS(c, "size_warning", sc_util_apen_get_size_warning, sc_util_apen_hold_size_warning);
    
    

	/* you CAN'T call this from the patcher */
	class_addmethod(c, (method)sc_util_apen_assist,			"assist",		A_CANT, 0);

	class_register(CLASS_BOX, c);
	sc_util_apen_class = c;
}

/////function definitions

void sc_util_apen_assist(t_sc_util_apen *x, void *b, long m, long a, char *s)
{
	if (m == ASSIST_INLET) { //inlet
        sprintf(s, "Inlet %ld: List of size %ld to add data to ApEn series / messages in", a, x->series_vector_size);
	}
	else {	// outlet
        if(a == 0) {
            sprintf(s, "Outlet %ld: Approximate Entropy", a);
        } else {
            sprintf(s, "Outlet %ld: dumpout", a);
        }
	}
}

void sc_util_apen_free(t_sc_util_apen *x)
{
    //start by freeing the test values
    if(x->test_value) {
        double* temp = x->test_value;
        for(int i = 0; i < x->series_max_length && temp; i++) {
            double* t2 = temp;
            temp++;
            sysmem_freeptr(t2);
        }
    }

}


t_max_err sc_util_apen_notify(t_sc_util_apen *x, t_symbol *s, t_symbol *msg, void *sender, void *data){
    t_symbol *attrname;
    
    if (msg == gensym("attr_modified")) {       // check notification type
        attrname = (t_symbol *)object_method((t_object *)data, gensym("getname"));      // ask attribute object for name
        object_post((t_object *)x, "changed attr name is %s",attrname->s_name);
    }
    return 0;
}

void sc_util_apen_getstate(t_sc_util_apen* x) {
    //output all attributes
    
    //pattern length
    void* state = sysmem_newptr(sizeof(t_atom) * 2);
    t_atom* pat_list = (t_atom*)state;
    t_atom* pat_temp = pat_list;
    atom_setsym(pat_temp, gensym("pattern_length"));
    pat_temp++;
    atom_setlong(pat_temp, x->pattern_length);
    outlet_list(x->out, gensym("pattern_length"), 2, (t_atom*)state);
    pat_temp = NULL;
    pat_list = NULL;
    
    //similarity
    t_atom* sim_list = (t_atom*)state;
    atom_setsym(sim_list, gensym("similarity"));
    sim_list++;
    atom_setfloat(sim_list, x->similarity);
    outlet_list(x->out, gensym("similarity"), 2, (t_atom*)state);
    sim_list = NULL;
    
    //calc on input
    t_atom* coi_list = (t_atom*)state;
    atom_setsym(coi_list, gensym("calculate_on_input"));
    coi_list++;
    atom_setlong(coi_list, x->calc_on_input);
    outlet_list(x->out, gensym("calculate_on_input"), 2, (t_atom*)state);
    coi_list = NULL;
    
    //max series length
    t_atom* temp_list = (t_atom*)state;
    atom_setsym(temp_list, gensym("series_length"));
    temp_list++;
    atom_setlong(temp_list, x->series_max_length);
    outlet_list(x->out, gensym("series_length"), 2, (t_atom*)state);
    
    //current series length
    temp_list = (t_atom*)state;
    atom_setsym(temp_list, gensym("current_size"));
    temp_list++;
    atom_setlong(temp_list, x->series_length);
    outlet_list(x->out, gensym("current_length"), 2, (t_atom*)state);
    
    //size warning
    temp_list = (t_atom*)state;
    atom_setsym(temp_list, gensym("size_warning"));
    temp_list++;
    atom_setlong(temp_list, x->hold_size_warning);
    outlet_list(x->out, gensym("size_warning"), 2, (t_atom*)state);
    
    //vector length
    temp_list = (t_atom*)state;
    atom_setsym(temp_list, gensym("vector_size"));
    temp_list++;
    atom_setlong(temp_list, x->series_vector_size);
    outlet_list(x->out, gensym("vector_size"), 2, (t_atom*)state);
    
    temp_list = NULL;
    sysmem_freeptr(state);
    
    sc_util_apen_dump(x);
    
}

void sc_util_apen_int(t_sc_util_apen *x, long n)
{
    
    critical_enter(0);
    
    if(x->series_length < x->series_max_length) {
        double* temp = x->test_value + x->series_length;
        *temp = (double) n;
        x->series_length++;
    } else {
        double* temp = x->test_value;
        double* temp2 = temp;
        temp2++;
        
        sysmem_copyptr(temp2, temp, sizeof(double) * (x->series_max_length - 1));
        temp = x->test_value + (x->series_max_length - 1);
        *temp = (double)n;

    }
    
    critical_exit(0);
    
    if(x->calc_on_input == 1) {
        sc_util_apen_calculate(x);
    }
}

void sc_util_apen_float(t_sc_util_apen *x, double f)
{
    critical_enter(0);
    
    if(x->series_length < x->series_max_length) {
        double* temp = x->test_value + x->series_length;
        
        *temp = f;
        x->series_length++;
    } else {
        double* temp = x->test_value;
        double* temp2 = temp;
        temp2++;
        
        sysmem_copyptr(temp2, temp, sizeof(double) * (x->series_max_length - 1));
        temp = x->test_value + (x->series_max_length - 1);
        *temp = f;
    }
    
    critical_exit(0);
    
    if(x->calc_on_input == 1) {
        sc_util_apen_calculate(x);
    }
}

void sc_util_apen_list(t_sc_util_apen *x, t_symbol* a, long argc, t_atom *argv) {
    

    t_atom* arg_temp = argv;
    long data_list_size = argc;
    long arg_offset = 0;
    if(argc > x->series_max_length)
    {
        data_list_size = x->series_max_length;
        arg_offset = argc - x->series_max_length;
    }
    
    double data_list[data_list_size];
    arg_temp += arg_offset;
    
    for(int i = 0; i < data_list_size && (i + arg_offset) < argc; i++, arg_temp++) {
        switch(atom_gettype(arg_temp)) {
            case A_LONG:
                data_list[i] = (double)atom_getlong(arg_temp);
                break;
            case A_FLOAT:
                data_list[i] = atom_getfloat(arg_temp);
                break;
            default:
                object_warn((t_object*)x, "Received non-numeric input");
                return;
        }
    }
    
    critical_tryenter(0);
    /*
    arg_temp = argv;
    int idx = 0;
    long data_size = argc;
    if(argc > x->series_max_length) {
        arg_temp += argc - x->series_max_length;
        idx = argc - x->series_max_length;
        data_size = x->series_max_length;
    }
    //double data_list;// = (double*)sysmem_newptr(sizeof(double) * data_size);
     
     for(; idx < argc; idx++, arg_temp++, data_temp++) {
        switch(atom_gettype(arg_temp)) {
            case A_LONG:
                *data_temp = (double)atom_getlong(arg_temp);
                break;
            case A_FLOAT:
                *data_temp = atom_getfloat(arg_temp);
                break;
        }
     }
    */
    
    //double* data_temp = data_list;

    
    long tot_size = x->series_length + data_list_size;
    
    long del_idx = 0;
    
    if(tot_size > x->series_max_length) {
        del_idx = tot_size - x->series_max_length;
        double* temp0 = x->test_value + del_idx;
        double* temp1 = x->test_value;
        sysmem_copyptr(temp0, temp1, sizeof(double) * (x->series_max_length - del_idx));
        x->series_length = (x->series_length - del_idx > 0) ? (x->series_length - del_idx) : 0;
    }
    
    double* temp = x->test_value + x->series_length;
    
    sysmem_copyptr(&data_list, temp, sizeof(double) * data_list_size);
    /*
    //free data
    data_temp = data_list;
    for(int i = 0; i < data_size; i++) {
        double* d2 = data_temp;
        data_temp++;
        sysmem_freeptr(d2);
    }
    */
    x->series_length += data_list_size;
    
    critical_exit(0);
    
    if(x->calc_on_input == 1) {
        sc_util_apen_calculate(x);
    }
}

void sc_util_apen_anything(t_sc_util_apen *x, t_symbol *s, long ac, t_atom *av)
{
    object_warn((t_object*)x, "Received something I'm unsure about");
}

//attempt calculating the current ApEn value without adding new data
void sc_util_apen_bang(t_sc_util_apen *x)
{
    sc_util_apen_calculate(x);
}


void sc_util_apen_dump(t_sc_util_apen *x) {
 
    critical_tryenter(0);

    if(x->series_length > 0){
        double* d = x->test_value;
        
        void* mem = sysmem_newptr(sizeof(t_atom) * (x->series_length + 1));
        t_atom* list = (t_atom*)mem;
        t_atom* temp_list = list;
        atom_setsym(temp_list, gensym("values"));
        temp_list++;
        for(int i = 0; i < x->series_length; i++, d++, temp_list++) {
            atom_setfloat(temp_list, *d);
        }
        outlet_list((void*)x->out, gensym("values"), x->series_length, list);
        
        sysmem_freeptr(mem);
        
    }
    critical_exit(0);
    
}

//empties list of data
void sc_util_apen_clear(t_sc_util_apen *x){
    
    critical_enter(0);

    double* temp = x->test_value;
    
    for(int i = 0; i < x->series_length; i++, temp++) {
        *temp = 0;
    }
    
    x->series_length = 0;
    
    critical_exit(0);
    
}

//sets the maximum length of the series
void sc_util_apen_set_series_length(t_sc_util_apen *x, void *attr, long argc, t_atom *argv){
    if(argc && argv) {
        
        long temp_sl = 0;
        
        switch(atom_gettype(argv)){
            case A_LONG:
                temp_sl = atom_getlong(argv);
                break;
            case A_FLOAT:
                temp_sl = (long)atom_getfloat(argv);
                break;
            default:
                object_error((t_object *)x, "bad value for series_length");
                return;
                break;
        }

        
        if(temp_sl > ((2 * x->pattern_length) + 1) && temp_sl != x->series_max_length) {
            critical_enter(0);
            
            double* temp = (double*)sysmem_newptr(sizeof(double) * temp_sl);
            long byte_length = 0;
            
            if(temp_sl > x->series_max_length) {
                byte_length = sizeof(double) * x->series_max_length;
            } else {
                byte_length = sizeof(double) * temp_sl;
            }
            sysmem_copyptr(x->test_value, temp, byte_length);
            
            if(x->series_length > temp_sl) {
                x->series_length = temp_sl;
            }
            
            x->series_max_length = temp_sl;
            
            critical_exit(0);
        } else if(temp_sl != x->series_max_length){
            object_error((t_object *)x, "Series length too short, must >= %d", (2 * x->pattern_length) + 1);
        }
    }
}

void sc_util_apen_get_series_length(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv){
    char alloc;
    long sl = 0;
    
    atom_alloc(argc, argv, &alloc);
    sl = x->series_max_length;
    atom_setlong(*argv, sl);
}

//sets the size of the data vector and clears the list
void sc_util_apen_set_vector_size(t_sc_util_apen *x, void *attr, long argc, t_atom *argv){
    if(argc && argv) {
        long temp_vs = 0;
        
        switch(atom_gettype(argv)){
            case A_LONG:
                temp_vs = atom_getlong(argv);
                break;
            case A_FLOAT:
                temp_vs = (double)atom_getfloat(argv);
                break;
            default:
                object_error((t_object *)x, "bad value received for vector_size");
                return;
                break;
        }
        if(temp_vs > 0) {
            x->series_vector_size = 1;
        } else {
            object_error((t_object *)x, "Vector Size must be a positive integer");
        }
    }
}

void sc_util_apen_get_vector_size(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv){
    char alloc;
    long vs = 0;
    
    atom_alloc(argc, argv, &alloc);
    vs = x->series_vector_size;
    atom_setlong(*argv, vs);
}

//sets the size of the pattern to be computed
void sc_util_apen_pattern_length(t_sc_util_apen *x, void *attr, long argc, t_atom *argv){
    if(argc && argv) {
        long temp_pl = 0;
        
        switch(atom_gettype(argv)){
            case A_LONG:
                temp_pl = atom_getlong(argv);
                break;
            case A_FLOAT:
                temp_pl = (double)atom_getfloat(argv);
                break;
            default:
                object_error((t_object *)x, "bad value received for pattern_length");
                return;
                break;
        }
        
        if(temp_pl <= (x->series_max_length / 2) - 1 && temp_pl > 1){
            x->pattern_length = temp_pl;
        } else if(temp_pl > (x->series_max_length / 2) - 1){
            object_error((t_object *)x, "pattern_length must be <= %d", (x->series_max_length / 2) - 1);
        } else {
            object_error((t_object *)x, "pattern_length must be an integer > 1");
        }
    }
}

void sc_util_apen_get_pattern_length(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv){
    char alloc;
    long pl = 0;
    
    atom_alloc(argc, argv, &alloc);
    pl = x->pattern_length;
    atom_setlong(*argv, pl);
}

//sets the threshold for pattern similarity
void sc_util_apen_similarity(t_sc_util_apen *x, void *attr, long argc, t_atom *argv){
    if(argv && argc) {
        double temp_sim = atom_getfloat(argv);
        
        if(temp_sim > 0.0) {
            x->similarity = temp_sim;
        } else {
            object_error((t_object *)x, "Similarity must be > 0.0, received %f", temp_sim);
        }
    }
}

void sc_util_apen_get_similarity(t_sc_util_apen *x,  t_object *attr, long *argc, t_atom **argv) {
    
    char alloc;
    double sim = 0.0;
    
    atom_alloc(argc, argv, &alloc);
    sim = x->similarity;
    atom_setfloat(*argv, sim);
}

//sets whether or not to attempt calculating ApEn when a new data point is received
void sc_util_apen_calc_on_input(t_sc_util_apen *x, void *attr, long argc, t_atom *argv){
    if(argv && argc) {
        long temp_coi = atom_getlong(argv);
        
        switch(atom_gettype(argv)) {
            case A_LONG:
                temp_coi = atom_getlong(argv);
                break;
            case A_FLOAT:
                temp_coi = (long)atom_getfloat(argv);
                break;
            default:
                object_error((t_object *)x, "bad value received for calculate_on_input");
                return;
                break;
        }
        
        if(temp_coi > 1){temp_coi = 1;}
        if(temp_coi < 0){temp_coi = 0;}
        x->calc_on_input = temp_coi;
    }
}

void sc_util_apen_get_coi(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv) {
    char alloc;
    long coi = 0;
    
    atom_alloc(argc, argv, &alloc);
    coi = x->calc_on_input;
    atom_setlong(*argv, coi);
}

void sc_util_apen_get_cur_size(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv){
    char alloc;
    long csize = 0;
    
    atom_alloc(argc, argv, &alloc);
    csize = x->series_length;
    atom_setlong(*argv, csize);
}

void sc_util_apen_set_cur_size(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv) {
    
}

void sc_util_apen_hold_size_warning(t_sc_util_apen *x, void *attr, long argc, t_atom *argv){
    if(argc && argv) {
        long temp_sw = 0;
        
        switch(atom_gettype(argv)) {
            case A_LONG:
                temp_sw = atom_getlong(argv);
                break;
            case A_FLOAT:
                temp_sw = (long)atom_getfloat(argv);
                break;
            default:
                object_error((t_object *)x, "bad value received for size_warning");
                return;
                break;
        }
        if(temp_sw >= 1) {temp_sw = 1;}
        if(temp_sw <= 0) {temp_sw = 0;}
        
        x->hold_size_warning = temp_sw;
    }
}


void sc_util_apen_get_size_warning(t_sc_util_apen *x, t_object *attr, long *argc, t_atom **argv){
    char alloc;
    long sw = 0;
    
    atom_alloc(argc, argv, &alloc);
    sw = x->hold_size_warning;
    atom_setlong(*argv, sw);
}


void *sc_util_apen_new(t_symbol *s, long argc, t_atom *argv)
{
	t_sc_util_apen *x = NULL;
    

	if ((x = (t_sc_util_apen *)object_alloc(sc_util_apen_class))) {
        //Set initial values
        x->calc_on_input = 1;
        x->hold_size_warning = 1;
        x->pattern_length = 3;
        x->series_length = 0;
        x->series_vector_size = 1; //for N-D vectors
        x->series_max_length = 50;
        x->similarity = 1.0;
		x->out = outlet_new(x, 0L);
        x->out2 = outlet_new(x, NULL);
        
        //allocate memory for the initial data series
        x->test_value = (double*)sysmem_newptr(sizeof(double) * x->series_max_length);
        
        //process arguments typed into object box
        attr_args_process(x, argc, argv);
        
    } else {
        poststring("Failed to create new ApEn");
    }
    
    
	return (x);
}

void sc_util_apen_calculate(t_sc_util_apen *x) {
    
    //check to make sure there is enough stored data to get meaningful results
    if(x->series_length < x->pattern_length * 2) {
        //check if the user has declined to have warnings sent to the console when there is insufficient data
        if(x->hold_size_warning == 1){ //warn user of insufficient data
            object_warn((t_object*)x, "Not enough data to calculate approximate entropy.");
            object_warn((t_object*)x, "Need %d data points, have %d", x->pattern_length * 2, x->series_length);
        }
        //exit function, do not attempt to calculate
        return;
    } else {
        
        //STEP 1 : Compute for pattern length
        double m0[x->series_length - x->pattern_length + 1]; //array for holding number of windows within similarity index of current window
        double avg_ratio0 = 0; //average number of windows within the similarity index for Cm(0...i)
        
        //temporary pointer to the data set
        double* temp = x->test_value;
        
        //outer loop for iterating through each possible window included in the data set of size m (pattern length)
        for(int i = 0; i < x->series_length - x->pattern_length + 1; i++, temp++) {
            m0[i] = 0.0; //initialize current index value
            double* temp2 = x->test_value; //second temporary pointer to data set
            //inner loop, iterate through all possible windows of size m to compare against the current window from the outer loop
            for(int j = 0; j < x->series_length - x->pattern_length + 1; j++, temp2++) {
                //compute the maximum distance between elements in both windows, add 1 to the index value if lees than or equal to similarity index
                m0[i] += (sc_util_apen_maxdist(temp, temp2, x->pattern_length, x->similarity) <= x->similarity) ? 1 : 0;
            }
            //get the percent of windows similar enough (total # of similar windows / total number of windows)
            m0[i] /= (x->series_length - x->pattern_length) + 1;
            //add to the average variable
            avg_ratio0 += m0[i];
        }
        //take the average percentage (Ci(m))
        avg_ratio0 /= (x->series_length - x->pattern_length) + 1;
        
        //STEP 2 : Compute for pattern length + 1
        double m1[x->series_length - (x->pattern_length + 1) + 1];
        double avg_ratio1 = 0;
        temp = x->test_value;
        for(int i = 0; i < x->series_length - (x->pattern_length + 1) + 1; i++, temp++) {
            m1[i] = 0;
            double* temp2 = x->test_value;
            for(int j = 0; j < x->series_length - (x->pattern_length + 1) + 1; j++, temp2++){
                m1[i] += (sc_util_apen_maxdist(temp, temp2, x->pattern_length + 1, x->similarity) <= x->similarity) ? 1 : 0;
            }
            m1[i] /= x->series_length - (x->pattern_length + 1) + 1;
            avg_ratio1 += m1[i];
        }
        
        //Ci(m+1)
        avg_ratio1 /= x->series_length - (x->pattern_length + 1) + 1;
        
        //STEP 3: get the Approximate Entropy as the Natural Logarithm of (Ci(m) / Ci(m)+1)
        double apen = log(avg_ratio0 / ((avg_ratio1 > 0.0) ? avg_ratio1 : 0.0000001)); //included a way to avoid division by 0 errors
        
        //outlet the value to the user
        outlet_float(x->out2, apen);
    }
}

//function for calculating the maximum pair-wise distance of members between two vectors
/* The function only compares members at matching indeces.
 Example:
 
 Vec1:          Vec2:           Dist:
 1     ->       5        =      4
 2     ->       4        =      2
 3     ->       3        =      0
 4     ->       2        =      2
 5     ->       1        =      4
 
 Maximum Distance : 4
 
 */

/* Paramters
 d0 - Double pointer to vector 0
 d1 - Double pointer to vector 1
 l - long size of both vectors, vectors should be the same length
 r - similarity index. used for optimizing when a distance is greater than the similarity index.
 */
double sc_util_apen_maxdist(double* d0, double* d1, long l, double r) {
    double md = 0.0;
    
    double* t = d0;
    double* t1 = d1;
    for(int i = 0; i < (l - 1); i++, t++, t1++) {
        double* t2 = t;
        t2++;
        for(int j = 1; j < (l - i); j++, t2++){
           md = (fabs(*t2 - *t) > md) ? md = fabs(*t2 - *t) : md;
            if(md > r) { //if the distance exceeds the similarity index r, just return the value, no need to calculate further
                return md;
            }
        }
    }
    
    return md;
}
