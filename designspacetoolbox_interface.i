%module dspace_interface

%{
#define SWIG_FILE_WITH_INIT
#include <DSStd.h>
#include <DSSSystem.h>
#include <DSTypes.h>
#include <DSDataSerialization.pb-c.h>
        
%}

/* Type Map Data */
%typemap(in) char ** {
        /* Check if is a list */
        if (PyList_Check($input)) {
                int size = PyList_Size($input);
                int i = 0;
                $1 = (char **) malloc((size+1)*sizeof(char *));
                for (i = 0; i < size; i++) {
                        PyObject *o = PyList_GetItem($input,i);
                        if (PyString_Check(o))
                        $1[i] = PyString_AsString(PyList_GetItem($input,i));
                        else {
                                PyErr_SetString(PyExc_TypeError,"list must contain strings");
                                free($1);
                                return NULL;
                        }
                }
                $1[i] = 0;
        } else {
                PyErr_SetString(PyExc_TypeError,"not a list");
                return NULL;
        }
}


%typemap(in) const DSUInteger * {
        /* Check if is a list */
        if (PyList_Check($input)) {
                int size = PyList_Size($input);
                int i = 0;
                $1 = (DSUInteger *) malloc((size+1)*sizeof(DSUInteger));
                for (i = 0; i < size; i++) {
                        PyObject *o = PyList_GetItem($input,i);
                        if (PyInt_Check(o))
                                $1[i] = PyInt_AsLong(PyList_GetItem($input,i));
                        else {
                                PyErr_SetString(PyExc_TypeError,"list must contain ints");
                                free($1);
                                return NULL;
                        }
                }
                $1[i] = 0;
        } else {
                PyErr_SetString(PyExc_TypeError,"not a list");
                return NULL;
        }
}

%typemap(out) char * {
        if ($1 == NULL) {
                return NULL;
        }
        $result = PyString_FromFormat("%s", $1);
        DSSecureFree($1);
}

%typemap(out) const char * {
        if ($1 == NULL) {
                return NULL;
        }
        $result = PyString_FromFormat("%s", $1);
}

%typemap(out) DSUInteger {
        $result = PyInt_FromLong((unsigned long)$1);
}

%typemap(in) DSUInteger {
        $1 = (DSUInteger) PyLong_AsUnsignedLongMask($input);
}

%typemap(in) const DSUInteger {
        $1 = (DSUInteger) PyLong_AsUnsignedLongMask($input);
}

%typemap(out) const DSUInteger {
        $result = PyInt_FromLong((unsigned long)$1);
}


%typemap(out) const DSVariable * {
        DSVariable * variable = NULL;
        variable = $1;
        PyObject * list = NULL;
        if (variable == NULL) {
                $result = NULL;
                return NULL;
        }
        list = PyList_New(2);
        PyList_SetItem(list, 0, PyString_FromFormat("%s", DSVariableName(variable)));
        PyList_SetItem(list, 1, PyFloat_FromDouble(DSVariableValue(variable)));
        $result = list;
}

%typemap(out) DSVertices *  {
        DSUInteger i, j;
        PyObject *tuple = NULL;
        DSVertices * vertices = $1;
        if (vertices == NULL) {
                $result = NULL;
                return NULL;
        }
        $result = PyList_New(vertices->numberOfVertices);
        for (i = 0; i < vertices->numberOfVertices; i++) {
                tuple = PyTuple_New(vertices->dimensions);
                for (j = 0; j < vertices->dimensions; j++) {
                        PyTuple_SetItem(tuple, j, PyFloat_FromDouble(DSVerticesVertexAtIndex(vertices, i)[j]));
                }
                PyList_SetItem($result, i, tuple);
        }
        DSVerticesFree(vertices);
}

%typemap(out) DSMatrix * {
        DSUInteger i, j;
        PyObject *tuple = NULL;
        DSMatrix *matrix = $1;
        if (matrix == NULL) {
                Py_RETURN_NONE;
//                $result = NULL;
//                return NULL;
        }
        $result = PyList_New(DSMatrixRows(matrix));
        for (i = 0; i < DSMatrixRows(matrix); i++) {
                tuple = PyTuple_New(DSMatrixColumns(matrix));
                for (j = 0; j < DSMatrixColumns(matrix); j++) {
                        PyTuple_SetItem(tuple, j, PyFloat_FromDouble(DSMatrixDoubleValue(matrix, i, j)));
                }
                PyList_SetItem($result, i, tuple);
        }
        DSMatrixFree(matrix);
}

%typemap(out) const DSMatrix * {
        DSUInteger i, j;
        PyObject *tuple = NULL;
        DSMatrix *matrix = $1;
        if (matrix == NULL) {
                Py_RETURN_NONE;
                //                $result = NULL;
                //                return NULL;
        }
        $result = PyList_New(DSMatrixRows(matrix));
        for (i = 0; i < DSMatrixRows(matrix); i++) {
                tuple = PyTuple_New(DSMatrixColumns(matrix));
                for (j = 0; j < DSMatrixColumns(matrix); j++) {
                        PyTuple_SetItem(tuple, j, PyFloat_FromDouble(DSMatrixDoubleValue(matrix, i, j)));
                }
                PyList_SetItem($result, i, tuple);
        }
//        DSMatrixFree(matrix);
}

%typemap(in) const DSCase ** {
        /* Check if is a list */
        if (PyList_Check($input)) {
                int size = PyList_Size($input);
                int i = 0;
                $1 = (const DSCase **) malloc((size+1)*sizeof(char *));
                for (i = 0; i < size; i++) {
                        PyObject *o = PyList_GetItem($input,i);
                        if (SWIG_ConvertPtr(o, &($1[i]), $descriptor(DSCase *), SWIG_POINTER_EXCEPTION) == -1) {
                                PyErr_SetString(PyExc_TypeError,"list must contain DSCase objects");
                                free($1);
                                return NULL;
                        }
                }
                $1[i] = 0;
        } else {
                PyErr_SetString(PyExc_TypeError,"not a list");
                return NULL;
        }
}

%include "/usr/local/include/designspace/DSStd.h"
%include "/usr/local/include/designspace/DSErrors.h"
%include "/usr/local/include/designspace/DSMemoryManager.h"
%include "/usr/local/include/designspace/DSVariable.h"
%include "/usr/local/include/designspace/DSMatrix.h"
%include "/usr/local/include/designspace/DSMatrixArray.h"
%include "/usr/local/include/designspace/DSExpression.h"
%include "/usr/local/include/designspace/DSGMASystem.h"
%include "/usr/local/include/designspace/DSSSystem.h"
%include "/usr/local/include/designspace/DSCase.h"
%include "/usr/local/include/designspace/DSDesignSpace.h"
%include "/usr/local/include/designspace/DSVertices.h"
%include "/usr/local/include/designspace/DSDictionary.h"
%include "/usr/local/include/designspace/DSStack.h"
%include "/usr/local/include/designspace/DSCyclicalCase.h"
%include "/usr/local/include/designspace/DSNVertexEnumeration.h"


%inline %{

typedef DSCase DSPseudoCase;
//extern DSCase * DSSWIGPseudoCaseAsCase(DSPseudoCase * pseudo)
//{
//        return pseudo;
//}
        
extern PyObject * DSSWIGDSCyclicalCaseEncodedBytes(DSCyclicalCase * aCase)
{
        PyObject * pyBuf = NULL;
        DSCyclicalCaseMessage * message;
        size_t length;
        unsigned char * buffer;
        message = DSCyclicalCaseEncode(aCase);
        if (message == NULL) {
                goto bail;
        }
        length = dscyclical_case_message__get_packed_size(message);
        buffer = DSSecureMalloc(sizeof(char)*length);
        dscyclical_case_message__pack(message, buffer);
        pyBuf = PyByteArray_FromStringAndSize(buffer, length);
        DSSecureFree(buffer);
        dscyclical_case_message__free_unpacked(message, NULL);
bail:
        return pyBuf;
}
        
extern DSCyclicalCase * DSSWIGDSCyclicalCaseDecodeFromByteArray(PyObject * byteArray)
{
        DSCyclicalCase * aCase = NULL;
        size_t length;
        const char * buffer;
        if (PyByteArray_Check(byteArray) == 0) {
                goto bail;
        }
        length = PyByteArray_Size(byteArray);
        buffer = PyByteArray_AsString(byteArray);
        aCase = DSCyclicalCaseDecode(length, buffer);
bail:
        return aCase;
}
        
extern PyObject * DSSWIGDSCaseEncodedBytes(DSCase * aCase)
{
        PyObject * pyBuf = NULL;
        DSCaseMessage * message;
        size_t length;
        unsigned char * buffer;
        message = DSCaseEncode(aCase);
        if (message == NULL) {
                goto bail;
        }
        length = dscase_message__get_packed_size(message);
        buffer = DSSecureMalloc(sizeof(char)*length);
        dscase_message__pack(message, buffer);
        pyBuf = PyByteArray_FromStringAndSize(buffer, length);
        DSSecureFree(buffer);
        dscase_message__free_unpacked(message, NULL);
bail:
        return pyBuf;
}
        
extern DSCase * DSSWIGDSCaseDecodeFromByteArray(PyObject * byteArray)
{
        DSCase * aCase = NULL;
        size_t length;
        const char * buffer;
        if (PyByteArray_Check(byteArray) == 0) {
                goto bail;
        }
        length = PyByteArray_Size(byteArray);
        buffer = PyByteArray_AsString(byteArray);
        aCase = DSCaseDecode(length, buffer);
bail:
        return aCase;
}
        
extern PyObject * DSSWIGDSDesignSpaceEncodedBytes(DSDesignSpace * ds) {
        PyObject * pyBuf = NULL;
        DSDesignSpaceMessage * message;
        size_t length;
        unsigned char * buffer;
        message = DSDesignSpaceEncode(ds);
        if (message == NULL) {
                goto bail;
        }
        length = dsdesign_space_message__get_packed_size(message);
        buffer = DSSecureMalloc(sizeof(char)*length);
        dsdesign_space_message__pack(message, buffer);
        pyBuf = PyByteArray_FromStringAndSize(buffer, length);
        DSSecureFree(buffer);
        dsdesign_space_message__free_unpacked(message, NULL);
bail:
        return pyBuf;
}

extern DSDesignSpace * DSSWIGDSDesignSpaceDecodeFromByteArray(PyObject * byteArray)
{
        DSDesignSpace * ds = NULL;
        size_t length;
        const char * buffer;
        if (PyByteArray_Check(byteArray) == 0) {
                goto bail;
        }
        length = PyByteArray_Size(byteArray);
        buffer = PyByteArray_AsString(byteArray);
        ds = DSDesignSpaceDecode(length, buffer);
bail:
        return ds;
}

extern DSDictionary * DSSWIGDSDictionaryFromPyDict(PyObject * pydict) {
        DSDictionary * dictionary = NULL;
        if (PyDict_Check(pydict) == false) {
                PyErr_SetString(PyExc_TypeError,"not a dictionary");
                goto bail;
        }
        int size = PyDict_Size(pydict);
        int i = 0;
        char * key_string, *value_string;
        PyObject * keys = PyDict_Keys(pydict);
        dictionary = DSDictionaryAlloc();
        for (i = 0; i < size; i++) {
                PyObject *key = PyList_GetItem(keys, i);
                PyObject *value = PyDict_GetItem(pydict, key);
                if (PyString_Check(key) && PyString_Check(value)) {
                        key_string = PyString_AsString(key);
                        value_string = PyString_AsString(value);
                        DSDictionaryAddValueWithName(dictionary, key_string, strdup(value_string));
                } else {
                        PyErr_SetString(PyExc_TypeError,"dictionary must contain strings");
                        DSDictionaryFreeWithFunction(dictionary, DSSecureFree);
                        dictionary = NULL;
                        break;
                }
        }
bail:
        return dictionary;
}
        
extern void DSSWIGDSDictionaryFreeCharValues(DSDictionary * dictionary)
{
        DSDictionaryFreeWithFunction(dictionary, DSSecureFree);
}
        
extern DSCyclicalCase * DSSWIGVoidAsSubcase(void *ptr)
{
        return ptr;
}

extern DSCase * DSSWIGVoidAsCase(void * ptr)
{
        return ptr;
}
        
extern DSVertices * DSSWIGVoidAsVertices(void * ptr)
{
        return ptr;
}
        
extern DSExpression * DSSWIGVoidAsExpression(void * ptr)
{
        return ptr;
}
        
extern DSDesignSpace * DSSWIGDesignSpaceParseWrapper(char ** const strings, const DSUInteger numberOfEquations, char ** Xd_list, const DSUInteger numberOfXd)
{
        DSUInteger i;
        DSVariablePool * Xd = DSVariablePoolAlloc();
        for (i = 0; i < numberOfXd; i++) {
                DSVariablePoolAddVariableWithName(Xd, Xd_list[i]);
        }
        DSDesignSpace * ds = DSDesignSpaceByParsingStrings(strings, Xd, numberOfEquations);
        DSVariablePoolFree(Xd);
        return ds;
}
        
extern DSDesignSpace * DSSWIGDesignSpaceParseWrapperWithXi(char ** const strings, const DSUInteger numberOfEquations, char ** Xd_list, const DSUInteger numberOfXd, char ** Xi_list, const DSUInteger numberOfXi)
{
        DSUInteger i;
        DSVariablePool * Xd = DSVariablePoolAlloc();
        DSVariablePool * Xi = DSVariablePoolAlloc();
        for (i = 0; i < numberOfXd; i++) {
                DSVariablePoolAddVariableWithName(Xd, Xd_list[i]);
        }
        for (i = 0; i < numberOfXi; i++) {
                DSVariablePoolAddVariableWithName(Xi, Xi_list[i]);
        }
        DSDesignSpace * ds = DSDesignSpaceByParsingStringsWithXi(strings, Xd, Xi, numberOfEquations);
        DSVariablePoolFree(Xd);
        DSVariablePoolFree(Xi);
        return ds;
}

extern DSGMASystem * DSSWIGGMASystemParseWrapper(char ** const strings, const DSUInteger numberOfEquations, char ** Xd_list, const DSUInteger numberOfXd)
{
        DSUInteger i;
        DSVariablePool * Xd = DSVariablePoolAlloc();
        for (i = 0; i < numberOfXd; i++) {
                DSVariablePoolAddVariableWithName(Xd, Xd_list[i]);
        }
        DSGMASystem * gma = DSGMASystemByParsingStrings(strings, Xd, numberOfEquations);
        DSVariablePoolFree(Xd);
        return gma;
}

extern DSSSystem * DSSWIGSSystemParseWrapper(char ** const strings, const DSUInteger numberOfEquations, char ** Xd_list, const DSUInteger numberOfXd)
{
        DSUInteger i;
        DSVariablePool * Xd = DSVariablePoolAlloc();
        for (i = 0; i < numberOfXd; i++) {
                DSVariablePoolAddVariableWithName(Xd, Xd_list[i]);
        }
        DSSSystem * ssys = DSSSystemByParsingStrings(strings, Xd, numberOfEquations);
        DSVariablePoolFree(Xd);
        return ssys;
}

extern DSExpression ** DSExpressionArrayFromVoid(void * pointer)
{
        return pointer;
}
        
extern DSExpression * DSExpressionAtIndexOfExpressionArray(DSExpression ** expressions, DSUInteger index)
{
        DSExpression * expression = NULL;
        expression = expressions[index];
        return expression;
}

extern DSUInteger DSUIntegerAtIndexOfIntegerArray(DSUInteger * array, DSUInteger index)
        {
                DSUInteger integer = NULL;
                integer = array[index];
                return integer;
        }

extern DSCase * DSCaseAtIndexOfArray(DSCase ** array, DSUInteger index)
{
        void * object = NULL;
        object = array[index];
        return object;
}
        
extern const char * DSDictionaryKeyAtIndex(const DSDictionary * dict, DSUInteger index)
{
        const char ** names = DSDictionaryNames(dict);
        DSUInteger i = 0;
        if (i >= DSDictionaryCount(dict)) {
                return NULL;
        }
        return names[index];
}

extern PyObject * DSSSystemPositiveRootsSWIG(const DSSSystem *ssys, const DSVariablePool *Xi0) {
        bool isMarginal = false;
        DSUInteger numberOfRoots;
        PyObject * list = PyList_New(2);
        numberOfRoots = DSSSystemPositiveRoots(ssys, Xi0, &isMarginal);
        PyList_SetItem(list, 0, PyInt_FromLong((long int)numberOfRoots));
        PyList_SetItem(list, 1, PyInt_FromLong((long int)isMarginal));
        return list;
}
        
%}


#define VERSION "0.01.1"