#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>


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
//  Function: EvaluateFieldAtLocation
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
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float EvaluateFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F)
{

    //I had to rewrite this entire equation because it kept breaking.
    int xIndexBefore, xIndexAfter, yIndexBefore, yIndexAfter;
    float xValueBefore, xValueAfter, yValueBefore, yValueAfter;
    for (int i=0; i<dims[0]; i++){
        if (pt[0] == X[i]){
            xIndexBefore = i;
            xIndexAfter = i;
            xValueBefore = X[i];
            xValueAfter = X[i];
            break;
        } else if ((pt[0] > X[i]) && (pt[0] < X[i+1])) {
            xIndexAfter = i+1;
            xIndexBefore = i;
            xValueBefore = X[i];
            xValueAfter = X[i+1];
            break;
        }
    }

    for (int j=0; j<dims[0]; j++){
        if (pt[1] == Y[j]){
            yIndexAfter = j;
            yIndexBefore = j;
            yValueBefore = Y[j];
            yValueAfter = Y[j];
            break;
        } else if ((pt[1] > Y[j]) && (pt[1] < Y[j+1])){
            yIndexBefore = j;
            yIndexAfter = j+1;
            yValueBefore = Y[j];
            yValueAfter = Y[j+1];
            break;
        }
    }

    float botLeftValue, botRightValue, topLeftValue, topRightValue;

    botLeftValue = F[(yIndexBefore * dims[0]) + xIndexBefore];
    botRightValue = F[(yIndexBefore * dims[0]) + xIndexAfter];
    topLeftValue = F[(yIndexAfter * dims[0]) + xIndexBefore];
    topRightValue = F[(yIndexAfter * dims[0]) + xIndexAfter];

    float xT, yT;
    float topInterpol, botInterpol, vertInterpol;

    if (topRightValue == topLeftValue){
        topInterpol = topLeftValue;
    } else {
        xT = (float)(pt[0] - xValueBefore)/(xValueAfter - xValueBefore);
        topInterpol = topLeftValue + (xT * (topRightValue - topLeftValue));
    }

    if (botRightValue == botLeftValue){
        botInterpol = botLeftValue;
    } else {
        xT = (float)(pt[0] - xValueBefore)/(xValueAfter - xValueBefore);
        botInterpol = botLeftValue + (xT * (botRightValue - botLeftValue));
    }

    if (topInterpol == botInterpol){
        vertInterpol = topInterpol;
    } else {
        yT = (float)(pt[1] - yValueBefore) / (yValueAfter - yValueBefore);
        vertInterpol = botInterpol + (yT * (topInterpol - botInterpol));
    }

    return vertInterpol;
}


void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
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

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue 
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,0.5) 
//        F=1: (1.0,1.0,1.0) 
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
//  Note: color map is in float (0 to 1) and output is in unsigned char
//        (0 to 255), so the floating point values you calculate for 
//        the three color channels need to be multiplied by 255 when
//        being assigned to RGB.
//
// ****************************************************************************

void ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
    unsigned char a, b;
    a = F * 255;
    b = ((F / 2) + .5) * 255;

    RGB[0] = a;
    RGB[1] = a;
    RGB[2] = b;
}


// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color 
//     using a difference colormap.
//
//     The difference color map has:
//        F=0: (0,0,0.5) 
//        F=0.5: (1.0,1.0,1.0) 
//        F=1: (0.5, 0, 0)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
//  Note: color map is in float (0 to 1) and output is in unsigned char
//        (0 to 255), so the floating point values you calculate for 
//        the three color channels need to be multiplied by 255 when
//        being assigned to RGB.
//
// ****************************************************************************

void ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
    // this fucking sucks as a method, but it works, so it's fine I guess
    unsigned char a, b, c;
    if (F <= 0.5){
        b = 255 * (2 * F);
        a = 255 * (2 * F);
        c = (0.5 + F) * 255;
    } else {
        a = 255 * (1.5 - F);
        b = 255 * (2.0 - (2 * F));
        c = 255 * (2.0 - (2 * F));
    }
    RGB[0] = a;
    RGB[1] = b;
    RGB[2] = c;
}

// ****************************************************************************
//  Function: ApplyHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0 <= F <= 1) to a color using 
//     an HSV rainbow colormap.
//
//     The rainbow colormap uses a saturation = 1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees.
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1  
//       RGB (output):  the location to store the color
//      
//  Note: as with the first two functions, make sure to multiple by 255 
//        when converting floats to unsigned chars.
//
// ****************************************************************************

void ApplyHSVColorMap(float F, unsigned char *RGB)  
{
    //I feel morally bad the quality of this code
    unsigned char a, b, c;
    float hue = F * 360;
    std::cout << "hue " << hue << "\n";
    if (hue <= 60) {
        a = 255;
        c = 0;
        b = 255 * (hue / 60);
        RGB[0] = a; 
        RGB[1] = b;
        RGB[2] = c;
        std::cout << " a " << +a << " b " << +b << " c " << +c << "\n"; 
        return;
    } else if ( hue <= 120) {
        b = 255;
        c = 0;
        a = 255 - (255 * ((hue - 60) / 60));
        RGB[0] = a; 
        RGB[1] = b;
        RGB[2] = c;
        std::cout << " a " << +a << " b " << +b << " c " << +c << "\n"; 
        return;
    } else if (hue <= 180) {
        b = 255;
        a = 0;
        c = 255 * ((hue - 120) / 60);
        RGB[0] = a; 
        RGB[1] = b;
        RGB[2] = c;
        std::cout << " a " << +a << " b " << +b << " c " << +c << "\n"; 
        return;
    } else if ( hue <= 240) {
        std::cout  << "hmmm " << ((hue - 180) / 60) * 255 << "\n";
        a = 0;
        c = 255;
        b = 255 - (((hue - 180) / 60) * 255);
        RGB[0] = a; 
        RGB[1] = b;
        RGB[2] = c;
        std::cout << " a " << +a << " b " << +b << " c " << +c << "\n";  
        return;
    } else if ( hue <= 300 ) {
        b = 0;
        c = 255;
        a = 255 * ((hue - 240) / 60);
        RGB[0] = a; 
        RGB[1] = b;
        RGB[2] = c;
        std::cout << " a " << +a << " b " << +b << " c " << +c << "\n"; 
        return;
    } else {
        a = 255;
        b = 0;
        c = 255 - (((hue - 300) / 60) * 255);
        RGB[0] = a; 
        RGB[1] = b;
        RGB[2] = c;
        std::cout << " a " << +a << " b " << +b << " c " << +c << "\n"; 
        return;
    }
}


int main()
{
    int  i, j;

    // creating a new reader object, naming the file to pull data from
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj3_data.vtk");
    rdr->Update();

    //pulling the dimensions of the data, rectilinear grid
    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    // vtk rectilinear grib object pulls the X coordinates
    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int nx = 500;
    int ny = 500;

    vtkImageData *images[3];
    unsigned char *buffer[3];
    for (i = 0 ; i < 3 ; i++)
    {
        images[i] = NewImage(nx, ny);
        buffer[i] = (unsigned char *) images[i]->GetScalarPointer(0,0,0);
    }

    for (i = 0 ; i < 3 ; i++)
       for (j = 0 ; j < 3*nx*ny ; j++)
            buffer[i][j] = 0;

    float xT, yT;
    float diff = 5.02 - 1.2;

    for (i = 0 ; i < nx ; i++){
        xT = (i /  499.0);
        for (j = 0 ; j < ny ; j++)
        {
            yT = j / 499.0;
            // ITERATE OVER PIXELS
            float pt[2];
            pt[0] = (18 * xT) - 9; //This is the interpolation formula, just rearranged a little. Also, I'm pretty sure it goes from -10 to 10;
            //riperooni, I should have not just tried to go off the dataset.
            pt[1] = (18 * yT) - 9;
            float f = EvaluateFieldAtLocation(pt, dims, X, Y, F);
            float normalizedF = (f - 1.2) / diff ; //...; see step 5 re 1.2->5.02
            
            // I TAKE OVER HERE0
            int offset = 3*(j*nx+i);
            ApplyBlueHotColorMap(normalizedF, buffer[0]+offset);
            ApplyDifferenceColorMap(normalizedF, buffer[1]+offset);
            ApplyHSVColorMap(normalizedF, buffer[2]+offset);
        }
    }

    WriteImage(images[0], "bluehot");
    WriteImage(images[1], "difference");
    WriteImage(images[2], "hsv");
}
