
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkDataSetWriter.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkUnsignedCharArray.h>


int main(int argc, char *argv[])
{
   vtkDataSetReader *reader = vtkDataSetReader::New();
   reader->SetFileName("noise.vtk");
   reader->Update();

   //create the plane
   vtkPlane *plane = vtkPlane::New();
   //create the cutter;
   vtkCutter *cutter = vtkCutter::New();
   //set the plane as the cut function
   cutter->SetCutFunction(plane);
   cutter->SetInputData(reader->GetOutput());

   vtkDataSetMapper *mapper = vtkDataSetMapper::New();
   mapper->SetInputConnection(cutter->GetOutputPort());
   
   vtkLookupTable *lut = vtkLookupTable::New();
   double *tmp = 0;
   vtkUnsignedCharArray *substitute = vtkUnsignedCharArray::New();
   //we iterate through, and create new color sets;
   substitute->SetNumberOfComponents(4);
   substitute->SetNumberOfTuples(256);
   for (int i=0; i<256; i++){
      substitute->InsertTuple4(i, int(255 * (255 - i) / 255), 0, int(255 * i / 255), 255);
      // this approach for interpolating came from the kitware examples in the vtk docs
   }
   //we reject their table, and substitute our own.
   lut->SetTable(substitute);
   mapper->SetLookupTable(lut);
   mapper->SetScalarRange(1,6);
   lut->Build();

   vtkActor *actor = vtkActor::New();
   actor->SetMapper(mapper);

   vtkRenderer *ren = vtkRenderer::New();
   ren->AddActor(actor);

   vtkRenderWindow *renwin = vtkRenderWindow::New();
   renwin->SetSize(768, 768);
   renwin->AddRenderer(ren);

   vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
   iren->SetRenderWindow(renwin);
   renwin->Render();
   iren->Start();
}


