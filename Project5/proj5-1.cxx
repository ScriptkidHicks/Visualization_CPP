/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}

// ****************************************************************************
//  Function: caseGenerate
//
//  Arguments:
//      xIndex: the x index of the bottom left corner
//      yIndex: the y index of the bottom right corner
//      X: an array of X values
//      Y: an array of Y values
//      F: an array of float values at points of the mesh
//      dims: the dimensions of the array
//      iso: the iso value. In this program it will be 3.2, but this function will have it passed in
//
//  Returns:  the int associated with which case we're in
//
// ****************************************************************************

int caseGenerate(int xIndex, int yIndex, float *X, float *Y, float *F, int *dims, float iso){
  int vert0, vert1, vert2, vert3;
  float tempFVal;

  int returnval = 0;
  int beforeind = (yIndex * dims[0]) + xIndex;

// ok, gotta figure out how to index into F correctly...
  tempFVal = F[beforeind];
  returnval = tempFVal >= iso ? 1 : 0;
  tempFVal = F[beforeind + dims[0]];
  returnval += tempFVal >= iso ? 2 : 0;
  tempFVal = F[beforeind + 1];
  returnval += tempFVal >= iso ? 4 : 0;
  tempFVal = F[beforeind + dims[0] + 1];
  returnval += tempFVal >= iso ? 8 : 0;

  // now that we have the 4 0 to 1 values, we can do some conversion into a single number;
  // actually, I think we can do this without weird binary work. Instead we use pseudobinary

  return returnval;
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}

// ****************************************************************************
//  Function: reverseInterpolate
//
//  Arguments:
//      a: the spacial value of location a
//      b: the spacial value of locaiton b
//      fOfA: the F field value at location A
//      fOfB: the F field value at location B
//      iso: the iso value, in this case 3.2, this is F(x)
//
//  Returns:  float ( the 'x' value which would produce the seen interpolation)
//
// ****************************************************************************

float reverseInterpolate(float a, float b, float fOfA, float fOfB, float iso){
  // To reverse intorpolate, just interpolate with reversed values.
  // there might be a faster way to do this? But my other method wasn't less
  // complicateed.
  float tmp1 = (iso - fOfA) / (fOfB - fOfA);
  float tmp2 = a + tmp1 * (b - a);
  return tmp2;
}


class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj5.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    SegmentList sl;
    sl.AddSegment(-10, -10, +10, -10); // Add segment (-10,-10) -> (+10, -10)
    sl.AddSegment(-10, +10, +10, +10);
    sl.AddSegment(-10, -10, -10, +10);
    sl.AddSegment(+10, -10, +10, +10);

// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!

    //ok, so I think we iterate over each cell here, and I think we have to normalize outputs on a size of -10 to 10
    
    // THINGS WE KNOW
    // the dimensions of both arrays run from -10 to +10, so there's no need to normalize.
    // I can just do a nested for loop, and not worry about the number of cells
    // we also know the isovalue is 3.2
    float isoValue = 3.2;
    int caseVal = 0;
    int beforeind;
    float edge0Point, edge1Point, edge2Point, edge3Point;

    for (int i=0; i<(dims[0] - 1); i++){
      for (int k=0; k<(dims[1] - 1); k++){
        //ok, here is where we figure out which case it is
        caseVal = caseGenerate(i, k, X, Y, F, dims, isoValue);
        beforeind = (k * dims[0]) + i;


        /*
        This version is considerably more verbose than the previous version I 
        turned in, at least in terms of code size, but there's a very good 
        reason for it. Though it might be less clear to look at, it means that 
        we are are not performing all 4 calculations very time through the loop,
        and only performing as many calculations as necessary. Previously I was 
        having the program calculate all 4 edge points even when it was cases 0
        or 15 (which can be a lot of wasted processing time when the vast 
        majority of cells are empty.)

        I'm working on a version that uses dictionary stored functions rather than a switch statement, meaning that the looking time for calling the appropriate function will be O(1)
        */

        switch (caseVal)
        {
        case 0: case 15:
          //do nothing;
          break;
        case 1: case 14:
          //calculate the endpoints for the relavant edges
          edge0Point = reverseInterpolate(X[i], X[i+1], F[beforeind], F[beforeind + 1], isoValue);
          edge1Point = reverseInterpolate(Y[k], Y[k+1], F[beforeind], F[beforeind + dims[0]], isoValue);
          // then place a line between those relevant edges
          sl.AddSegment(X[i], edge1Point, edge0Point, Y[k]);
          break;

        case 2: case 13:
          edge1Point = reverseInterpolate(Y[k], Y[k+1], F[beforeind], F[beforeind + dims[0]], isoValue);
          edge2Point = reverseInterpolate(X[i], X[i+1], F[beforeind + dims[0]], F[beforeind + dims[0] + 1], isoValue);
          sl.AddSegment(X[i], edge1Point, edge2Point, Y[k+1]);
          break;

        case 3: case 12:
          edge0Point = reverseInterpolate(X[i], X[i+1], F[beforeind], F[beforeind + 1], isoValue);
          edge2Point = reverseInterpolate(X[i], X[i+1], F[beforeind + dims[0]], F[beforeind + dims[0] + 1], isoValue);
          sl.AddSegment(edge0Point, Y[k], edge2Point, Y[k+1]);
          break;
        case 4: case 11:
          edge0Point = reverseInterpolate(X[i], X[i+1], F[beforeind], F[beforeind + 1], isoValue);
          edge3Point = reverseInterpolate(Y[k], Y[k+1], F[beforeind + 1], F[beforeind + 1 + dims[0]], isoValue);
          sl.AddSegment(edge0Point, Y[k], X[i+1], edge3Point);
          break;
        case 5: case 10:
          edge1Point = reverseInterpolate(Y[k], Y[k+1], F[beforeind], F[beforeind + dims[0]], isoValue);
          edge3Point = reverseInterpolate(Y[k], Y[k+1], F[beforeind + 1], F[beforeind + 1 + dims[0]], isoValue);
          sl.AddSegment(X[i], edge1Point, X[i+1], edge3Point);
          break;
        case 6:
          edge0Point = reverseInterpolate(X[i], X[i+1], F[beforeind], F[beforeind + 1], isoValue);
          edge1Point = reverseInterpolate(Y[k], Y[k+1], F[beforeind], F[beforeind + dims[0]], isoValue);
          edge2Point = reverseInterpolate(X[i], X[i+1], F[beforeind + dims[0]], F[beforeind + dims[0] + 1], isoValue);
          edge3Point = reverseInterpolate(Y[k], Y[k+1], F[beforeind + 1], F[beforeind + 1 + dims[0]], isoValue);
          sl.AddSegment(edge2Point, Y[k+1], X[i+1], edge3Point);
          sl.AddSegment(X[i], edge1Point, edge0Point, Y[k]);
          break;
        case 7: case 8:
          edge2Point = reverseInterpolate(X[i], X[i+1], F[beforeind + dims[0]], F[beforeind + dims[0] + 1], isoValue);
          edge3Point = reverseInterpolate(Y[k], Y[k+1], F[beforeind + 1], F[beforeind + 1 + dims[0]], isoValue);
          sl.AddSegment(edge2Point, Y[k+1], X[i+1], edge3Point);
          break;
        case 9:
          edge0Point = reverseInterpolate(X[i], X[i+1], F[beforeind], F[beforeind + 1], isoValue);
          edge1Point = reverseInterpolate(Y[k], Y[k+1], F[beforeind], F[beforeind + dims[0]], isoValue);
          edge2Point = reverseInterpolate(X[i], X[i+1], F[beforeind + dims[0]], F[beforeind + dims[0] + 1], isoValue);
          edge3Point = reverseInterpolate(Y[k], Y[k+1], F[beforeind + 1], F[beforeind + 1 + dims[0]], isoValue);
          sl.AddSegment(X[i], edge1Point, edge2Point, Y[k+1]);
          sl.AddSegment(edge0Point, Y[k], X[i+1], edge3Point);
          break;
        default:
          break;
        }
      }
    }

    


    vtkPolyData *pd = sl.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
