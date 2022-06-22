
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

   //create the contour filter
   vtkContourFilter *filter = vtkContourFilter::New();
   filter->SetNumberOfContours(2);
   filter->SetValue(1, 2.4);
   filter->SetValue(2, 4);
   filter->SetInputConnection(reader->GetOutputPort());

   //write the data to a buffer
   vtkDataSetWriter *writer = vtkDataSetWriter::New();
   writer->SetFileName("contouredOutput.vtk");
   writer->SetInputConnection(filter->GetOutputPort());
   writer->Write();

   // ok, at this point we're going to need a new reader
   vtkDataSetReader *myreader = vtkDataSetReader::New();
   myreader->SetFileName("contouredOutput.vtk");
   myreader->Update();

   //now we create our own mapper
   vtkDataSetMapper *mymapper = vtkDataSetMapper::New();
   mymapper->SetInputData(myreader->GetOutput());

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
   mymapper->SetLookupTable(lut);
   mymapper->SetScalarRange(1,6);
   lut->Build();

   vtkActor *actor = vtkActor::New();
   actor->SetMapper(mapper);

   //create our own actor
   vtkActor *myactor = vtkActor::New();
   myactor->SetMapper(mymapper);

   vtkRenderer *ren = vtkRenderer::New();
   ren->AddActor(actor);
   
   //this time we create our own renderer
   vtkRenderer *myren = vtkRenderer::New();
   myren->AddActor(myactor);
   myren->SetViewport(0.5, 0, 1, 1);
   ren->SetViewport(0, 0, 0.5, 1);

   vtkRenderWindow *renwin = vtkRenderWindow::New();
   renwin->SetSize(1536, 768);
   renwin->AddRenderer(ren);
   renwin->AddRenderer(myren);

   vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
   iren->SetRenderWindow(renwin);
   renwin->Render();
   iren->Start();
}


