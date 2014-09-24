%module dspace_interface

%{
#define SWIG_FILE_WITH_INIT
#include <DSStd.h>
#include <DSSSystem.h>
#include <DSTypes.h>
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

/* Type Map Data */

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

//
//
///**
// * DSVariablePool functions available to the internal python module.
// */
//extern DSUInteger DSVariablePoolNumberOfVariables(const DSVariablePool *pool);
//extern DSVariablePool * DSVariablePoolAlloc(void);
//extern DSVariablePool * DSVariablePoolCopy(const DSVariablePool * const pool);
//extern void DSVariablePoolFree(DSVariablePool *pool);
//extern void DSVariablePoolAddVariableWithName(DSVariablePool *pool, const char * name);
//extern void DSVariablePoolSetValueForVariableWithName(const DSVariablePool *pool, const char *name, const double value);
//extern bool DSVariablePoolHasVariableWithName(const DSVariablePool *pool, const char * const name);
//extern double DSVariablePoolValueForVariableWithName(const DSVariablePool *pool, const char *const name);
//extern DSUInteger DSVariablePoolIndexOfVariableWithName(const DSVariablePool *pool, const char *name);
//extern void DSVariablePoolPrint(const DSVariablePool * const pool);
//extern DSMatrix * DSVariablePoolValuesAsVector(const DSVariablePool *pool, const bool rowVector);
//
//extern const DSVariable * DSVariablePoolVariableAtIndex(const DSVariablePool *pool, const DSUInteger index);
//extern void DSVariablePoolSetReadWrite(DSVariablePool *pool);
//extern void DSVariablePoolSetReadWriteAdd(DSVariablePool *pool);
//
//
///**
// * DSDesignSpace functions available to internal python module.
// */
//extern DSDesignSpace * DSDesignSpaceByParsingStrings(char * const * const strings, const DSVariablePool * const Xd_a, const DSUInteger numberOfEquations);
//void DSDesignSpaceFree(DSDesignSpace * ds);
//
//extern const DSDictionary * DSDesignSpaceSubcaseDictionary(const DSDesignSpace *ds);
//
//extern void DSDesignSpacePrint(const DSDesignSpace * ds);
//extern void DSDesignSpaceCalculateValidityOfCases(DSDesignSpace *ds);
//
//extern const DSVariablePool * DSDesignSpaceXi(const DSDesignSpace *ds);
//
//extern const DSUInteger DSDesignSpaceNumberOfEquations(const DSDesignSpace *ds);
//extern DSExpression ** DSDesignSpaceEquations(const DSDesignSpace *ds);
//extern const DSUInteger * DSDesignSpaceSignature(const DSDesignSpace *ds);
//
//extern const DSUInteger DSDesignSpaceNumberOfValidCases(const DSDesignSpace *ds);
//extern const DSUInteger DSDesignSpaceNumberOfCases(const DSDesignSpace *ds);
//
//extern DSCase * DSDesignSpaceCaseWithCaseNumber(const DSDesignSpace * ds, const DSUInteger caseNumber);
//extern DSCase * DSDesignSpaceCaseWithCaseSignature(const DSDesignSpace * ds, const DSUInteger * signature);
//
//extern const bool DSDesignSpaceCaseWithCaseNumberIsValid(const DSDesignSpace *ds, const DSUInteger caseNumber);
//extern const bool DSDesignSpaceCaseWithCaseSignatureIsValid(const DSDesignSpace *ds, const DSUInteger * signature);
//
////extern const DSStack * DSDesignSpaceSubcasesForCaseNumber(DSDesignSpace *ds, const DSUInteger caseNumber);
//extern const DSGMASystem * DSDesignSpaceGMASystem(const DSDesignSpace * ds);
//
//extern DSCase ** DSDesignSpaceCalculateAllValidCases(DSDesignSpace *ds);
//extern DSDictionary * DSDesignSpaceCalculateAllValidCasesForSlice(DSDesignSpace *ds, const DSVariablePool *lower, const DSVariablePool *upper);
//
//extern void DSDesignSpaceCalculateCyclicalCases(DSDesignSpace *ds);
//
///**
// * DSGMASystem functions available to internal python module
// */
//
//extern DSGMASystem * DSGMASystemByParsingStrings(char * const * const strings, const DSVariablePool * const Xd_a, const DSUInteger numberOfEquations);
//extern void DSGMASystemFree(DSGMASystem * gma);
//
//extern const DSUInteger DSGMASystemNumberOfEquations(const DSGMASystem *gma);
//extern DSExpression ** DSGMASystemEquations(const DSGMASystem *gma);
//
//extern const DSVariablePool *DSGMASystemXd(const DSGMASystem *gma);
//extern const DSVariablePool *DSGMASystemXd_a(const DSGMASystem *gma);
//extern const DSVariablePool *DSGMASystemXi(const DSGMASystem *gma);
//
//extern void DSGMASystemPrintEquations(const DSGMASystem *gma);
//
///**
// * DSSSystem functions available to internal python module
// */
//
//extern DSSSystem * DSSSystemByParsingStringList(char * const * const string, const DSVariablePool * const Xd_a, ...);
//extern DSSSystem * DSSSystemByParsingStrings(char * const * const strings, const DSVariablePool * const Xd_a, const DSUInteger numberOfEquations);
//
//extern DSSSystem * DSSSystemByRemovingAlgebraicConstraints(const DSSSystem * originalSSystem);
//extern double DSSSystemLogarithmicGain(const DSSSystem *ssys, const char *XdName, const char *XiName);
//
//
//extern const DSUInteger DSSSystemNumberOfEquations(const DSSSystem *ssys);
//extern DSExpression ** DSSSystemEquations(const DSSSystem *ssys);
//extern DSExpression ** DSSSystemSolution(const DSSSystem *ssys);
//extern DSExpression ** DSSSystemLogarithmicSolution(const DSSSystem *ssys);
//
//extern double DSSSystemSteadyStateFunction(const DSSSystem *ssys, const DSVariablePool *Xi0, const char * function);
//extern DSMatrix * DSSSystemSteadyStateValues(const DSSSystem *ssys, const DSVariablePool *Xi0);
//extern DSMatrix * DSSSystemSteadyStateFluxForDependentVariables(const DSSSystem * ssys,
//                                                                const DSVariablePool * Xd0,
//                                                                const DSVariablePool * Xi0);
//extern DSMatrix * DSSSystemSteadyStateFlux(const DSSSystem *ssys, const DSVariablePool *Xi0);
//extern DSMatrix * DSSSystemRouthArrayForPoolTurnover(const DSSSystem *ssys, const DSMatrix * F);
//extern DSMatrix * DSSSystemRouthArrayForSteadyState(const DSSSystem *ssys,
//                                                    const DSVariablePool *Xd0,
//                                                    const DSVariablePool *Xi0);
//extern DSMatrix * DSSSystemRouthArray(const DSSSystem *ssys, const DSVariablePool *Xi0);
//extern DSUInteger DSSSystemNumberOfPositiveRootsForRouthArray(const DSMatrix *routhArray);
//extern DSUInteger DSSSystemPositiveRootsForSteadyState(const DSSSystem *ssys,
//                                                       const DSVariablePool *Xd0,
//                                                       const DSVariablePool *Xi0);
//extern DSUInteger DSSSystemPositiveRoots(const DSSSystem *ssys, const DSVariablePool *Xi0);
//extern DSUInteger DSSSystemRouthIndex(const DSSSystem *ssys, const DSVariablePool *Xi0);
//extern DSUInteger DSSSystemCharacteristicEquationCoefficientIndex(const DSSSystem *ssys, const DSVariablePool *Xi0);
//
//
//extern void DSSSystemPrint(DSSSystem * ssys);
//extern void DSSSystemPrintEquations(DSSSystem * ssys);
//extern void DSSSystemPrintSolution(DSSSystem * ssys);
//extern void DSSSystemPrintLogarithmicSolution(DSSSystem *ssys);
//
//extern void DSSSystemFree(DSSSystem * ssys);
//
//extern const DSVariablePool * DSSSystemXd(const DSSSystem * const ssys);
//extern const DSVariablePool * DSSSystemXd_a(const DSSSystem * const ssys);
//extern const DSVariablePool * DSSSystemXi(const DSSSystem * const ssys);
//
///**
// * DSCase functions available to internal python module.
// */
//
//extern void DSCaseFree(DSCase *aCase);
//extern void DSCasePrint(const DSCase *aCase);
//extern DSExpression ** DSCaseEquations(const DSCase *aCase);
//
//extern DSExpression ** DSCaseSolution(const DSCase *aCase);
//extern DSExpression ** DSCaseLogarithmicSolution(const DSCase *aCase);
//
//extern DSExpression ** DSCaseConditions(const DSCase *aCase);
//extern DSExpression ** DSCaseLogarithmicConditions(const DSCase *aCase);
//
//extern DSExpression ** DSCaseBoundaries(const DSCase *aCase);
//extern DSExpression ** DSCaseLogarithmicBoundaries(const DSCase *aCase);
//extern DSUInteger DSCaseNumberOfEquations(const DSCase *aCase);
//extern DSUInteger DSCaseNumberOfConditions(const DSCase *aCase);
//
//extern DSVertices * DSCaseVerticesFor2DSlice(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const char * xVariable, const char *yVariable);
//extern DSVertices * DSCaseVerticesFor1DSlice(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const char * xVariable);
//extern DSVertices * DSCaseBoundingRangeForVariable(const DSCase *aCase, const char * variable);
//extern DSVertices * DSCaseBoundingRangeForVariableWithConstraints(const DSCase *aCase, const char * variable, DSVariablePool * lowerBounds, DSVariablePool * upperBounds);
//
//extern DSUInteger DSCaseNumber(const DSCase * aCase);
//extern const DSUInteger * DSCaseSignature(const DSCase * aCase);
//extern char * DSCaseSignatureToString(const DSCase *aCase);
//
//extern const bool DSCaseIsValid(const DSCase *aCase);
//extern const bool DSCaseIsValidAtSlice(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds);
//
//extern const DSSSystem * DSCaseSSystem(const DSCase *aCase);
//
//extern const bool DSCaseIntersectionIsValid(const DSUInteger numberOfCases, const DSCase **cases);
//extern const bool DSCaseIntersectionIsValidAtSlice(const DSUInteger numberOfCases, const DSCase **cases,  const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds);
//extern const bool DSCaseIsValidAtPoint(const DSCase *aCase, const DSVariablePool * variablesToFix);
//
//extern DSVertices * DSCaseIntersectionVerticesForSlice(const DSUInteger numberOfCases, const DSCase **cases, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const DSUInteger numberOfVariables, const char ** variables);
//
//extern const bool DSCaseIntersectionExceptSliceIsValid(const DSUInteger numberOfCases, const DSCase **cases, const DSUInteger numberOfExceptions, const char ** exceptionVarNames);
//extern const bool DSCaseIntersectionExceptSliceIsValidAtSlice(const DSUInteger numberOfCases, const DSCase **cases, const DSUInteger numberOfExceptions, const char ** exceptionVarNames, const DSVariablePool * lowerBounds, const DSVariablePool * upperBounds);
//extern DSVariablePool * DSCaseIntersectionExceptSliceValidParameterSet(const DSUInteger numberOfCases, const DSCase **cases, const DSUInteger numberOfExceptions, const char ** exceptionVarNames);
//extern DSVariablePool * DSCaseIntersectionExceptSliceValidParameterSetAtSlice(const DSUInteger numberOfCases, const DSCase **cases, const DSUInteger numberOfExceptions, const char ** exceptionVarNames, const DSVariablePool * lowerBounds, const DSVariablePool * upperBounds);
//
//extern double DSCaseLogarithmicGain(const DSCase *aCase, const char *XdName, const char *XiName);
//
//extern DSVariablePool * DSCaseValidParameterSet(const DSCase *aCase);
//extern DSVariablePool * DSCaseValidParameterSetAtSlice(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds);
///**
//* Subcases
//*/
//extern DSCyclicalCase * DSCyclicalCaseForCaseInDesignSpace(const DSDesignSpace * ds, const DSCase * aCase);
//extern void DSCyclicalCaseFree(DSCyclicalCase * aSubcase);
//
//extern const bool DSCyclicalCaseIsValid(const DSCyclicalCase *aSubcase);
//extern const bool DSCyclicalCaseIsValidAtSlice(const DSCyclicalCase *aSubcase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds);
//
//extern const DSDesignSpace * DSCyclicalCaseInternalDesignSpace(const DSCyclicalCase * subcase);
//extern const DSCase * DSCyclicalCaseOriginalCase(const DSCyclicalCase * cyclicalCase);
//
//extern DSCase * DSCyclicalCaseSubcaseWithCaseNumber(const DSCyclicalCase * cyclicalCase, const DSUInteger subcaseNumber);
//
//extern DSDictionary * DSCyclicalCaseCalculateAllValidSubcasesForSlice(const DSCyclicalCase * cyclicalCase,
//                                                                      const DSVariablePool *lower,
//                                                                      const DSVariablePool *upper);
//
//extern DSDictionary * DSCyclicalCaseVerticesForSlice(const DSCyclicalCase *cyclicalCase,
//                                                     const DSVariablePool * lowerBounds,
//                                                     const DSVariablePool *upperBounds,
//                                                     const DSUInteger numberOfVariables,
//                                                     const char ** variables);
//
//extern DSDictionary * DSCyclicalCaseVerticesFor2DSlice(const DSCyclicalCase *cyclicalCase,
//                                                       const DSVariablePool * lowerBounds,
//                                                       const DSVariablePool *upperBounds,
//                                                       const char * xVariable,
//                                                       const char *yVariable);
///**
// * Dictionary Functions
// */
//
//extern DSDictionary * DSDictionaryAlloc(void);
//extern void DSDictionaryFree(DSDictionary * aDictionary);
//extern void DSDictionaryFreeWithFunction(DSDictionary * aDictionary, void * freeFunction);
//
//extern DSUInteger DSDictionaryCount(const DSDictionary *aDictionary);
//extern void * DSDictionaryValueForName(const DSDictionary *dictionary, const char *name);
//
//extern const char ** DSDictionaryNames(const DSDictionary *aDictionary);
//
//extern void DSDictionaryAddValueWithName(DSDictionary *dictionary, const char * name, void *value);
//
//
///**
// * Vertices 
// */
////extern void DSVerticesFree(DSVertices *vertices);
//
////extern const bool DSVerticesAddVertex(DSVertices *vertices, const double * coordinates);
//
////extern const bool DSVerticesAreEqual(const DSVertices *vert1, const DSVertices *vert2);
//
////extern const double * DSVerticesVertexAtIndex(const DSVertices *vertices, const DSUInteger index);
//
//#if defined(__APPLE__) && defined(__MACH__)
//#pragma mark - Utility functions
//#endif
//
//
//extern void DSSecureFree(void * ptr);
//extern DSDictionary * DSDictionaryFromArray(void * array, DSUInteger size);
//extern DSExpression * DSExpressionByCompressingConstantVariables(const DSExpression *expression, const DSVariablePool * assumedConstant);
//extern DSExpression * DSExpressionByParsingString(const char *string);
//extern double DSExpressionEvaluateWithVariablePool(const DSExpression *expression, const DSVariablePool *pool);
//extern char * DSExpressionAsString(const DSExpression *expression);
//extern void DSExpressionFree(DSExpression *expression);
//
////extern DSCase * DSSWIGVoidAsCase(void * ptr);
//
////extern void * DSSWIGDesignSpaceParseWrapper(const DSVariablePool * const Xd, char ** const strings, const DSUInteger numberOfEquations);

%inline %{

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