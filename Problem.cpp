#include "Problem.hpp"
#include <math.h>
#include <cmath>

/*
		i = 0, j = Jmax	      i = Imax, j = Jmax
		+--------------------
		|
		|
		|
		|
		|
		+---------------------

		i = 0, j=0			  i = Imax, j = 0
*/


Problem::Problem(unsigned int N, unsigned int M):
_N(N), // Number of internal nodes
_M(M),
_Imax(N + 1), // Max index over i <-> csi
_Jmax(M + 1),
_deltaCsi(1.0/(_Imax)), // eta, csi space span [0,1]x[0,1] space
_deltaEta(1.0/(_Jmax)),
_deltaX(_L/_Imax),
_deltaY(_H/_Jmax),
_x(_Imax + 1, _Jmax + 1),
_y(_Imax + 1, _Jmax + 1),
_alpha(_Imax - 1, _Jmax - 1, 0.0),
_beta (_Imax - 1, _Jmax - 1, 0.0),
_gamma(_Imax - 1, _Jmax - 1, 0.0),
_d    (_Imax - 1, _Jmax - 1, 0.0),
_a1   (_Imax - 1, _Jmax - 1, 0.0),
_a2   (_Imax - 1, _Jmax - 1, 0.0),
_a3   (_Imax - 1, _Jmax - 1, 0.0),
_a4   (_Imax - 1, _Jmax - 1, 0.0),
_a5   (_Imax - 1, _Jmax - 1, 0.0),
_a6   (_Imax - 1, _Jmax - 1, 0.0),
_a7   (_Imax - 1, _Jmax - 1, 0.0),
_a8   (_Imax - 1, _Jmax - 1, 0.0)
{
	std::cout << "Problem parameters " << std::endl;
	std::cout << "deltaEta = " << _deltaEta << std::endl;
	std::cout << "deltaCsi = " << _deltaCsi << std::endl;
	std::cout << "Number of internal nodes along csi : " << N << std::endl;
	std::cout << "Number of internal nodes along eta : " << M << std::endl;
}

Problem::~Problem()
{
}

double pow2(double f)
{
	return f*f;
}

// Initialize with no bump
// void Problem::initialize()
// {
//     // Simple uniform rectangle in physical space
//     for (unsigned int i = 0; i <= _Imax; ++i)
//     {
//         double xI = i * _deltaX;  // with _deltaX = L / _Imax

//         for (unsigned int j = 0; j <= _Jmax; ++j)
//         {
//             double yJ = j * _deltaY;  // with _deltaY = H / _Jmax

//             _x(i,j) = xI;
//             _y(i,j) = yJ;
//         }
//     }
// }


void Problem::initialize()
{
	// x and y are set to zero everywhere except on the boundary 
	// Adjust boundary conditions for j = _Jmax ---> x = _L
	// and for j = 0 ---> x = 0.
	// Initialize internal points with algebraic grid

	std::cout << "Initializing x " << std::endl;

	// Set value of x at the left
	_x.setRow(0, 0);

	bool leadingEdge = true;
	bool trailingEdge = true;
	
	// Along csi (i)
	for (unsigned int i = 1; i < _Imax ; ++i)
	{
		// Step from one node to the other
		double xI = i * _deltaX;

		// bottom
		_x(i,0)     = _x(i - 1, 0) + _deltaX;
		// top
		_x(i,_Jmax) = _x(i - 1,_Jmax) + _deltaX;
				
	}
	
	// Fill interior
	for (unsigned int i = 1; i < _Imax; i++)
	{
		/* code */
		for (unsigned int j = 1; j < _Jmax ; ++j)
		{
			//_x(i,j) = 0.0 ;
			_x(i,j) = _x(i,j-1);
		
		}
	}
	
	_x.setRow(_Imax, _L);

	std::cout << "Initializing y " << std::endl;

	// Apply boundary conditions at j = 0 (bottom) and j = _Jmax (top)
	// If csi is between 2 and 3, use sinusoidal profile,
	// otherwise bottom = 0, top = _H

	for (unsigned int i = 0; i <= _Imax; ++i)
	{
		double xI = i * _deltaX;

	//	if (1<0)
		if (xI > 2.0 && xI < 3.0)
		{
			// Use x(i,0) / x(i,_Jmax) just like you used x(0,j) before
			_y(i,0)      = 0.1 * std::sin(M_PI * (_x(i,0)      - 2.0));       // bottom bump
			
			_y(i,_Jmax)  = 1.0 - 0.1 * std::sin(M_PI * (_x(i,_Jmax) - 2.0)); // top bump
			// _y(i,_Jmax)  = _H;
		}
		else
		{
			_y(i,0)      = 0.0;   // flat bottom
			_y(i,_Jmax)  = _H;    // flat top
		}

		// _y(i,0)      = 0.0;   // flat bottom
		// _y(i,_Jmax)  = _H;    // flat top
	}


	// Fill interior nodes in j by linear interpolation between bottom and top
	for (unsigned int i = 0; i <= _Imax; ++i)
	{
		for (unsigned int j = 1; j < _Jmax; ++j) // j = 1.._Jmax-1
		{
			double s = static_cast<double>(j) / static_cast<double>(_Jmax); // 0..1
			_y(i,j) = (1.0 - s) * _y(i,0) + s * _y(i,_Jmax);
		}
	}
}

void Problem::solve()
{
    unsigned int iter = 0;
    Matrix<double> xPrev(_Imax+1, _Jmax+1, 0.0);
    Matrix<double> yPrev(_Imax+1, _Jmax+1, 0.0);

    double errX, errY;

    do
    {
        // 1) update coefficients for this iteration
        updateCoeff();

        double errSumX = 0.0;
        double errSumY = 0.0;

        // 2) Gauss–Seidel for x on interior points
        for (unsigned int i = 1; i < _Imax; ++i)
        {
            for (unsigned int j = 1; j < _Jmax; ++j)
            {
                unsigned int ii = i - 1;
                unsigned int jj = j - 1;

                double xOld = _x(i,j);

                _x(i,j) =
                      _a1(ii,jj) * _x(i+1,j)
                    + _a2(ii,jj) * _x(i-1,j)
                    + _a3(ii,jj) * _x(i+1,j+1)
                    + _a4(ii,jj) * _x(i-1,j+1)
                    + _a5(ii,jj) * _x(i+1,j-1)
                    + _a6(ii,jj) * _x(i-1,j-1)
                    + _a7(ii,jj) * _x(i,  j+1)
                    + _a8(ii,jj) * _x(i,  j-1);

                errSumX += std::abs(_x(i,j) - xOld);
            }
        }

        // 3) Gauss–Seidel for y on interior points
        for (unsigned int i = 1; i < _Imax; ++i)
        {
            for (unsigned int j = 1; j < _Jmax; ++j)
            {
                unsigned int ii = i - 1;
                unsigned int jj = j - 1;

                double yOld = _y(i,j);

                _y(i,j) =
                      _a1(ii,jj) * _y(i+1,j)
                    + _a2(ii,jj) * _y(i-1,j)
                    + _a3(ii,jj) * _y(i+1,j+1)
                    + _a4(ii,jj) * _y(i-1,j+1)
                    + _a5(ii,jj) * _y(i+1,j-1)
                    + _a6(ii,jj) * _y(i-1,j-1)
                    + _a7(ii,jj) * _y(i,  j+1)
                    + _a8(ii,jj) * _y(i,  j-1);

                errSumY += std::abs(_y(i,j) - yOld);
            }
        }

		enforceBC(BOTTOM);
		enforceBC(TOP);

        errX = errSumX;
        errY = errSumY;

        std::cout << "-------------------------------------\n";
        std::cout << "iter = " << iter  << "\n";
        std::cout << "sum error for x : " << errX << "\n";
        std::cout << "sum error for y : " << errY << "\n";
        std::cout << "-------------------------------------\n";

        ++iter;
        if (iter > 100000)
        {
            std::cout << "Not converged !!!!!!!!!!!!!!!!!!!!!!\n";
            break;
        }

    } while (errX > _tol || errY > _tol);
}



void Problem::updateCoeff()
{
    const double eps = 0.0;

    for (unsigned int i = 1; i < _Imax; ++i)      // i = 1..Imax-1
    {
        for (unsigned int j = 1; j < _Jmax; ++j)  // j = 1..Jmax-1
        {	
			// alpha beta and gamma are defined only in internal points
            unsigned int ii = i - 1;
            unsigned int jj = j - 1;

            // alpha, beta, gamma
            double alpha = 1.0/(4.0 * pow2(_deltaEta)) *
                ( pow2(_x(i, j+1) - _x(i, j-1)) +
                  pow2(_y(i, j+1) - _y(i, j-1)) ) + eps;

            double beta  = 1.0/(4.0 * _deltaEta * _deltaCsi) *
                ( (_x(i+1, j) - _x(i-1, j)) * (_x(i, j+1) - _x(i, j-1)) +
                  (_y(i+1, j) - _y(i-1, j)) * (_y(i, j+1) - _y(i, j-1)) ) + eps;

            double gamma = 1.0/(4.0 * pow2(_deltaCsi)) *
                ( pow2(_x(i+1, j) - _x(i-1, j)) +
                  pow2(_y(i+1, j) - _y(i-1, j)) ) + eps;

            double d = 2.0*alpha/pow2(_deltaCsi) + 2.0*gamma/pow2(_deltaEta);

            _alpha(ii, jj) = alpha;
            _beta (ii, jj) = beta;
            _gamma(ii, jj) = gamma;
            _d    (ii, jj) = d;

            _a1(ii, jj) =  alpha / (pow2(_deltaEta) * d);
            _a2(ii, jj) =  alpha / (pow2(_deltaEta) * d);
            _a3(ii, jj) = -beta  / (2.0 * _deltaEta * _deltaCsi * d);
            _a4(ii, jj) =  beta  / (2.0 * _deltaEta * _deltaCsi * d);
            _a5(ii, jj) =  beta  / (2.0 * _deltaEta * _deltaCsi * d);
            _a6(ii, jj) = -beta  / (2.0 * _deltaEta * _deltaCsi * d);
            _a7(ii, jj) =  gamma / (pow2(_deltaCsi) * d);
            _a8(ii, jj) =  gamma / (pow2(_deltaCsi) * d);
        }
    }
}

void Problem::enforceBC(enum topBottom side)
{
    // Pick boundary and interior rows based on which side
    unsigned int jBoundary = 0;
    unsigned int jInterior = 1;

    if (side == TOP)
    {
        jBoundary = _Jmax;
        jInterior = _Jmax - 1;
    }
    else // BOTTOM
    {
        jBoundary = 0;
        jInterior = 1;
    }

    // index for tracking closest points to leading/trailing edge
    unsigned int idxX2 = 0;
    unsigned int idxX3 = _Imax;

    // Initialize with current distances  TODO: generalize leading and trailing edge
    double bestDist2 = std::abs(_x(idxX2, jBoundary) - 2.0);
    double bestDist3 = std::abs(_x(idxX3, jBoundary) - 3.0);

    // Loop over interior lines (horizontal)
    for (unsigned int i = 1; i < _Imax; ++i)
    {
        // slope on interior horizontal line using central difference
        double dx = _x(i+1, jInterior) - _x(i-1, jInterior);
        double dy = _y(i+1, jInterior) - _y(i-1, jInterior);

        // Avoid division by zero
        if (std::abs(dx) < 1e-14)
            continue;

        double slope = dy / dx;

        // Neumann orthogonality is imposed using boundary and interior
        double x1 = _x(i, jBoundary);
        double y1 = _y(i, jBoundary);
        double x2 = _x(i, jInterior);
        double y2 = _y(i, jInterior);

        // skip flat walls
        if (std::abs(slope) > 1e-14)
        {
            // Intersection between:
            //   y = y1 + slope * (x - x1)           (tangent through boundary point)
            //   y = y2 - (1/slope) * (x - x2)       (normal through interior point)
            double s  = slope;
            double is = 1.0 / s;

            double xNew = (x1 * s + x2 * is + y2 - y1) / (s + is);

            _x(i, jBoundary) = xNew;
        }

        // Copy the new x to snap the surface
        double xB = _x(i, jBoundary);

        // Apply Dirichlet to determine y on this boundary
        if (xB < 2.0 || xB > 3.0)
        {
            if (side == BOTTOM)
            {
                _y(i, jBoundary) = 0.0;
            }
            else // TOP
            {
                _y(i, jBoundary) = _H;
            }
        }
        else
        {
            if (side == BOTTOM)
            {
                _y(i, jBoundary) = 0.1 * std::sin(M_PI * (xB - 2.0));
            }
            else // TOP: bump into domain from y = _H
            {
                _y(i, jBoundary) = _H - 0.1 * std::sin(M_PI * (xB - 2.0));
            }
        }

        // Distance from leading edge
        double d2 = std::abs(xB - 2.0);
        if (d2 < bestDist2)
        {
            bestDist2 = d2;
            idxX2 = i;
        }

        // Distance from trailing edge
        double d3 = std::abs(xB - 3.0);
        if (d3 < bestDist3)
        {
            bestDist3 = d3;
            idxX3 = i;
        }
    }

    // Fix boundary points at exact leading / trailing edge
    _x(idxX2, jBoundary) = 2.0;
    _x(idxX3, jBoundary) = 3.0;

    if (side == BOTTOM)
    {
        _y(idxX2, jBoundary) = 0.0;
        _y(idxX3, jBoundary) = 0.0;
    }
    else // TOP
    {
        _y(idxX2, jBoundary) = _H;
        _y(idxX3, jBoundary) = _H;
    }

    // Optionally fix left corner exactly (applies to both top and bottom)
    _x(0, jBoundary) = 0.0;
    if (side == BOTTOM)
    {
        _y(0, jBoundary) = 0.0;
    }
    else
    {
        _y(0, jBoundary) = _H;
    }
}

void Problem::save()
{
	_y.save("y");
	_x.save("x");

	writeVTK(_x, _y, "output.vtk");
}

