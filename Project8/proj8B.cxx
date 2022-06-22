
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkDataSetWriter.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkLookupTable.h>


int main(int argc, char *argv[])
{
   vtkDataSetReader *reader = vtkDataSetReader::New();
   reader->SetFileName("noise.vtk");
   reader->Update();

   //contour the data
   vtkContourFilter *cf = vtkContourFilter::New();
   cf->SetNumberOfContours(2);
   cf->SetValue(1, 2.4);
   cf->SetValue(2, 4);
   cf->SetInputConnection(reader->GetOutputPort());

   //now we write the data to a buffer
   vtkDataSetWriter *writer = vtkDataSetWriter::New();
   writer->SetFileName("contouredOutput.vtk");
   writer->SetInputConnection(cf->GetOutputPort());
   writer->Write();

   //change what the reader points at
   reader->SetFileName("contouredOutput.vtk");
   reader->Update();

   vtkDataSetMapper *mapper = vtkDataSetMapper::New();
   mapper->SetInputData(reader->GetOutput());
   
   vtkLookupTable *lut = vtkLookupTable::New();
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

   cf->Delete();
   writer->Delete();
   reader->Delete();
   renwin->Delete();
   actor->Delete();
   ren->Delete();
   mapper->Delete();
   lut->Delete();
   iren->Delete();
}


