#ifndef vtkPVDataRepresentation_h
#define vtkPVDataRepresentation_h

#include <vtkUnstructuredGridAlgorithm.h>
#include "vtkIOCCMModule.h"

class VTKIOCCM_EXPORT vtkCCMReader : public vtkUnstructuredGridAlgorithm
{
public:
    static vtkCCMReader* New();
    vtkTypeMacro(vtkCCMReader, vtkUnstructuredGridAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent) override;
public:
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

protected:
    vtkCCMReader();
    ~vtkCCMReader();
    int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
        vtkInformationVector* outputVector) override;

    char* FileName;
private:
    vtkCCMReader(const vtkCCMReader&) = delete;
    void operator=(const vtkCCMReader&) = delete;

};


#endif