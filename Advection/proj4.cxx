#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkDataSetWriter.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSphereSource.h>

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
//  Function: EvaluateVectorFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//        The first value in the field is the x-component for the first point.
//        The second value in the field is the y-component for the first point.
//
//     rv (output): the interpolated field value. (0,0) if the location is out of bounds.
//
// ****************************************************************************

void EvaluateVectorFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F, float *rv)
{
    // IMPLEMENT ME!

    int xLessIndex, xMoreIndex, yLessIndex, yMoreIndex;
    float xLessValue, xMoreValue, yLessValue, yMoreValue;

    for (int i=0; i<dims[0]; i++){
        if (X[i] == pt[0]){
            xLessIndex = i;
            xMoreIndex = i;
            xLessValue = X[i];
            xMoreValue = X[i];
            break;
        } else if ((pt[0] > X[i]) && (pt[0] < X[i + 1])){
            xLessIndex = i;
            xMoreIndex = i+1;
            xLessValue = X[i];
            xMoreValue = X[i+1];
            break;
        }
    }

    for (int i=0; i<dims[1]; i++){
        if (Y[i] == pt[1]){
            yLessIndex = i;
            yMoreIndex = i;
            yLessValue = Y[i];
            yMoreValue = Y[i];
            break;
        } else if ((pt[1] > Y[i]) && (pt[1] < Y[i + 1])){
            yLessIndex = i;
            yMoreIndex = i+1;
            yLessValue = Y[i];
            yMoreValue = Y[i+1];
            break;
        }
    }

    //ok, we should now have the left and right index / mesh location values
    //Now we procure the values for each corner, but we're also going to need to procure the X/Y values (naming is gonna get confusing here)

    float botLeftXVal, botLeftYVal, topLeftXVal, topLeftYVal, topRightXVal, topRightYVal, botRightXVal, botRightYVal;

    botLeftXVal = F[2 * (yLessIndex * dims[0] +  xLessIndex)];
    botLeftYVal = F[2 * (yLessIndex * dims[0] + xLessIndex) + 1];
    //std::cout << "bot left y " << botLeftYVal << endl;

    botRightXVal = F[2 * (yLessIndex * dims[0] + xMoreIndex)];
    botRightYVal = F[2 * (yLessIndex * dims[0] + xMoreIndex) + 1];
    //std::cout << "bot right y " << botRightYVal << endl;

    topLeftXVal = F[2 * (yMoreIndex * dims[0] + xLessIndex)];
    topLeftYVal = F[2 * (yMoreIndex * dims[0] + xLessIndex) + 1];
    //std::cout << "top left y " << topLeftYVal << endl;

    topRightXVal = F[2 * (yMoreIndex * dims[0] + xMoreIndex)];
    topRightYVal = F[2 * (yMoreIndex * dims[0] + xMoreIndex)  + 1];
    //std::cout << "top right y " << topRightYVal << endl;

    float xT, yT;

    float topXInterpol, topYInterpol, botXInterpol, botYInterpol, vertXInterpol, vertYInterpol;

    if (topRightXVal == topLeftXVal){
        topXInterpol = topLeftXVal;
        topYInterpol = topLeftYVal;
    } else {
        xT = (float)(pt[0] - xLessValue)/(xMoreValue - xLessValue);
        topXInterpol = topLeftXVal + (xT * (topRightXVal - topLeftXVal));
        topYInterpol = topLeftYVal + (xT * (topRightYVal - topLeftYVal));
        botXInterpol = botLeftXVal + (xT * (botRightXVal - botLeftXVal));
        botYInterpol = botLeftYVal + (xT * (botRightYVal - botLeftYVal));
    }

    if (topLeftXVal == botLeftXVal){
        vertXInterpol = topLeftXVal;
        vertYInterpol = topLeftYVal;
    } else {
        yT = (float)(pt[1] - yLessValue)/(yMoreValue - yLessValue);
        vertXInterpol = botXInterpol + (yT * (topXInterpol - botXInterpol));
        vertYInterpol = botYInterpol + (yT * (topYInterpol - botYInterpol));
    }

    rv[0] = vertXInterpol; // setting the x-component of the velocity
    rv[1] = vertYInterpol; // setting the y-component of the velocity

}

// ****************************************************************************
//  Function: AdvectWithEulerStep
//
//  Arguments:
//     pt: the seed location (two-dimensions)
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//     h: The size of the Euler step
//     nsteps: The number of Euler steps to take
//     output_locations (output): An array of size 2*(nsteps+1).  It's first entry
//        should be the seed location.  The second entry should be the result
//        of the first advection step.  The final entry should be the result
//        of the final advection step.
//
// ****************************************************************************

void
AdvectWithEulerStep(const float *pt, const int *dims, const float *X, 
                    const float *Y, const float *F, 
                    float h, int nsteps, float *output_locations)
{
    // I am still convinced that I could do this in fewer lines, but I'm trying to figure out how
    float placement[2] = {0, 0};
    float combo[2] = {0, 0};
    output_locations[0] = pt[0];
    output_locations[1] = pt[1];
    for (int i=1; i<=nsteps; i++){
        placement[0] = output_locations[(i * 2) - 2];
        placement[1] = output_locations[(i * 2) - 1];
        EvaluateVectorFieldAtLocation(placement, dims, X, Y, F, combo);
        output_locations[i * 2] = combo[0] * h + placement[0];
        output_locations[(i * 2) + 1] = combo[1] * h + placement[1];
    }

}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// VTK files are only 3D, so the vector data is all of the form (X,Y,0).
// Remove the 0's since it is counter-intuitive for students who are 
// thinking of this as 2D data.
float *
Convert3DVectorDataTo2DVectorData(const int *dims, const float *F)
{
    float *rv = new float[dims[0]*dims[1]*2];
    int index3D = 0;
    int index2D = 0;
    for (int i = 0 ; i < dims[0] ; i++)
       for (int j = 0 ; j < dims[1] ; j++)
       {
           rv[index2D]   = F[index3D];
           rv[index2D+1] = F[index3D+1];
           index2D += 2;
           index3D += 3;
       }

    return rv;
}

/*
This is the RK4 version of the same feature Euler Advection implements.

*/
void AdvectWithRK4Step( const float *pt, const int *dims, const float *X, 
                        const float *Y, const float *F, 
                        float h, int nsteps, float *output_locations){
    output_locations[0] = pt[0];
    output_locations[1] = pt[1];

    float halfH = h/2;
    float sixthH = h/6;
    float k1[2] = {0, 0};
    float k2[2] = {0, 0};
    float k3[2] = {0, 0};
    float k4[2] = {0, 0};
    float combo[2] = {0, 0};
    float initialPlacement[2] = {0, 0};
    float k1Placement[2] = {0, 0};
    float k2Placement[2] = {0, 0};
    float k3Placement[2] = {0, 0};
    for (int i=1; i<=nsteps; i++){
        initialPlacement[0] = output_locations[(i * 2) - 2];
        initialPlacement[1] = output_locations[(i * 2) - 1];
        //now we calculalte k1;
        EvaluateVectorFieldAtLocation(initialPlacement, dims, X, Y, F, k1); 
        k1Placement[0] = initialPlacement[0] + ((halfH) * k1[0]);
        k1Placement[1] = initialPlacement[1] + ((halfH) * k1[1]);

        EvaluateVectorFieldAtLocation(k1Placement, dims, X, Y, F, k2);
        k2Placement[0] = initialPlacement[0] + ((halfH) * k2[0]);
        k2Placement[1] = initialPlacement[1] + ((halfH) * k2[1]);
        EvaluateVectorFieldAtLocation(k2Placement, dims, X, Y, F, k3);
        std::cout << "k3: " << k3[0] << " " << k3[1] << endl; 
        k3Placement[0] = initialPlacement[0] + (h * k3[0]);
        k3Placement[1] = initialPlacement[1] + (h * k3[1]);
        EvaluateVectorFieldAtLocation(k3Placement, dims, X, Y, F, k4); 
        output_locations[i * 2] = initialPlacement[0] + ((sixthH) * (k1[0] + (2 * k2[0]) + (2 * k3[0]) + k4[0]));
        output_locations[(i * 2) + 1] = initialPlacement[1] + ((sixthH) * (k1[1] + (2 * k2[1]) + (2 * k3[1]) + k4[1]));
    }

}

void PrintSteps(const char *solver, int nsteps, float *locations)
{
   cerr << "Printing output for solver " << solver << endl;
   for (int j = 0 ; j < nsteps+1 ; j++)
   {
       cerr << j << ": (" << locations[2*j] << ", " << locations[2*j+1] << ")" << endl;
   }
}

int main()
{
    int  i, j;

    // HANK'S CODE TO SET THINGS UP -- DON'T MODIFY THIS
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj4_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    if (dims[0] <= 0 || dims[1] <= 0)
    {
        cerr << "Was not able to successfully open file \"proj4_data.vtk\"" << endl;
        exit(EXIT_FAILURE);
    }
    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F_3D = (float *) rgrid->GetPointData()->GetVectors()->GetVoidPointer(0);
    float *F = Convert3DVectorDataTo2DVectorData(dims, F_3D);
    
    float seed[2] = { 1, -5 };
    
    // SANITY CHECK TO MAKE SURE VECTOR FIELD EVALUATION IS WORKING
    float vec[2];
    EvaluateVectorFieldAtLocation(seed, dims, X, Y, F, vec);
    cerr << "Velocity at (" << seed[0] <<", " << seed[1] << ") is (" << vec[0] << ", " << vec[1] << ")" << endl;

    float h = 0.01;
    const int nsteps = 100;
    float *output_locations = new float[2*(nsteps+1)];
    AdvectWithEulerStep(seed, dims, X, Y, F, h, nsteps, output_locations);
    PrintSteps("Euler", nsteps, output_locations);

    // Uncomment this code if you want to do the RK4 part of the project
    float *RK4_output_locations = new float[2*(nsteps+1)];
    AdvectWithRK4Step(seed, dims, X, Y, F, h, nsteps, RK4_output_locations);
    PrintSteps("RK4", nsteps, RK4_output_locations);
}
