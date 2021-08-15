/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * Flow around a 2D cylinder inside a channel, with the creation of a von
 * Karman vortex street. This example makes use of bounce-back nodes to
 * describe the shape of the cylinder. The outlet is modeled through a
 * Neumann (zero velocity-gradient) condition.
 */

#include "palabos2D.h"
#include "palabos2D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
///ʹ��D2Q9ģ��
#define DESCRIPTOR D2Q9Descriptor

/// Velocity on the parabolic Poiseuille profile ��������ٶȱ������� poiseuilleVelocity;
/// IncomprFlowParam �࣡���²���ѹ����������ֵ������ͨ����������IncomprFlowParam����������
/// plint��plb�ռ�Palabosר�Ŷ���ģ�������int
T poiseuilleVelocity(plint iY, IncomprFlowParam<T> const& parameters) {
//    T y = (T)iY / parameters.getResolution();
 //   T y = (T)iY /(parameters.getResolution() * parameters.getLy());
//    return 4.*parameters.getLatticeU() * (y-y*y);
    return parameters.getLatticeU();
}

/// Linearly decreasing pressure profile ѹ��
T poiseuillePressure(plint iX, IncomprFlowParam<T> const& parameters) {
    T Lx = parameters.getNx()-1;
    T Ly = parameters.getNy()-1;
    return 8.*parameters.getLatticeNu()*parameters.getLatticeU() / (Ly*Ly) * (Lx/(T)2-(T)iX);
}

/// Convert pressure to density according to ideal gas law ѹ��ת��Ϊ�����ܶ�
T poiseuilleDensity(plint iX, IncomprFlowParam<T> const& parameters) {
    return poiseuillePressure(iX,parameters)*DESCRIPTOR<T>::invCs2 + (T)1;
}

/// A functional, used to initialize the velocity for the boundary conditions
/// Array ���ݸ�ֵ
///ͨ����������࣬��poiseuilleVelocity��IncromprFlowParam�������
///����ʹ��operator�����ǻ�������virtual bool operator

template<typename T>
class PoiseuilleVelocity {
public:
//PoiseuilleVelocity��Ĺ��캯�����������������ȫ��ͬ�ģ�һ������ĳ�Ա������������ÿ�δ�������¶���ʱִ�С�
//���캯��Ҳ���Դ��в����������ڴ�������ʱ�ͻ�����󸳳�ʼֵ
    PoiseuilleVelocity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
	/// operatorֻ�����ٶ�
    void operator()(plint iX, plint iY, Array<T,2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

/// A functional, used to initialize a pressure boundary to constant density
template<typename T>
class ConstantDensity {
public:
    ConstantDensity(T density_)
        : density(density_)
    { }
    T operator()(plint iX, plint iY) const {
        return density;
    }
private:
    T density;
};

/// A functional, used to create an initial condition for the density and velocity��ʼ������
template<typename T>
class PoiseuilleVelocityAndDensity {
public:
    PoiseuilleVelocityAndDensity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, T& rho, Array<T,2>& u) const {
        rho = poiseuilleDensity(iX,parameters);
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

/// A functional, used to instantiate bounce-back nodes at the locations of the cylinder
/// Բ������ȫ���ص��߽磡
/// virtual��̬���麯��
template<typename T>
class CylinderShapeDomain2D : public plb::DomainFunctional2D {
public:
    CylinderShapeDomain2D(plb::plint cx_, plb::plint cy_, plb::plint radius)
        : cx(cx_),
          cy(cy_),
          radiusSqr(plb::util::sqr(radius))
    { }
    virtual bool operator() (plb::plint iX, plb::plint iY) const {
        return plb::util::sqr(iX-cx) + plb::util::sqr(iY-cy) <= radiusSqr;
    }
    virtual CylinderShapeDomain2D<T>* clone() const {
        return new CylinderShapeDomain2D<T>(*this);
    }
private:
    plb::plint cx;
    plb::plint cy;
    plb::plint radiusSqr;
};


void defineCylinderGeometry( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                             IncomprFlowParam<T> const& parameters,
                             OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition,
                             Array<plint,2> forceIds )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
	
	Box2D inlet(0, 0, 0, ny-1);
    Box2D outlet(nx-1,nx-1, 0, ny-1);

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet, boundary::outflow );
	
    setBoundaryVelocity ( lattice, inlet, PoiseuilleVelocity<T>(parameters) );
    setBoundaryDensity (
            lattice, outlet,
            ConstantDensity<T>(1.) );
    initializeAtEquilibrium (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocityAndDensity<T>(parameters) );
	
	
	//***ԭģ�����ûص��߽�
    // defineDynamics(lattice, lattice.getBoundingBox(),
                   // new CylinderShapeDomain2D<T>(cx,cy,radius),
                   // new plb::BounceBack<T,DESCRIPTOR>);
	//*******
/*/ Instead of plain BounceBack, use the dynamics MomentumExchangeBounceBack,
    //   to compute the momentum exchange, and thus, the drag and the lift on
    //   the obstacle, locally during collision.
    defineDynamics(lattice, lattice.getBoundingBox(),
                   new CylinderShapeDomain2D<T>(cx,cy,radius),
                   new MomentumExchangeBounceBack<T,DESCRIPTOR>(forceIds));
    initializeMomentumExchange(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(cx, cy, radius) );
*/
	MultiScalarField2D<bool> boolMask(parameters.getNx(), parameters.getNy());

    plb_ifstream ifile("geometry.dat");
    ifile >> boolMask;
	defineDynamics(lattice,  boolMask, new MomentumExchangeBounceBack<T,DESCRIPTOR>(forceIds), true);
	initializeMomentumExchange(lattice, lattice.getBoundingBox());
//  ������
    lattice.initialize();
}

void writeGif(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                               *computeVelocityNorm(lattice) );
}

void writeVTK(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", dx/dt);
}


int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
///  ����plbInit����ͷǿ�Ƶ��ã���֤˳��Ͳ��г�����һ�£�
    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 20.,  // Re
            40,       // N  �ֱ��ʣ����ߴ�Ϊ1��101������
            6.,        // lx
            1.         // l
    );
    const T logT     = (T)0.02;
//    const T imSave   = (T)0.06;
    const T vtkSave  = (T)1.;
    const T maxT     = (T)200.1;
	const T Tolerror = (T)1e-10;
	T old_Energy = (T)1e-3;

    writeLogFile(parameters, " flow");
	
///  �����New�������������Ͷ�̬�����ڴ棡��Ӧdelete�ͷ��ڴ棬��Ϊ�޷���ǰԤ֪��Ҫ�����ڴ����洢ĳ����������е��ض���Ϣ�������ڴ�Ĵ�С��Ҫ������ʱ����ȷ��
///  �������
	MultiBlockLattice2D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(),
            new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );
lattice.initialize();
    // OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        // boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    // cylinderSetup(lattice, parameters, *boundaryCondition);
// The drag and lift acting on the obstacle are computed with help of the
    //   internal statistics object of the lattice. For this purpose, they
    //   need to be registered first, as it is done in the following lines.
    Array<plint,2> forceIds;
	lattice.toggleInternalStatistics(true);
    forceIds[0] = lattice.internalStatSubscription().subscribeSum();
    forceIds[1] = lattice.internalStatSubscription().subscribeSum();

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createInterpBoundaryCondition2D<T,DESCRIPTOR>();
	
	
    defineCylinderGeometry(lattice, parameters, *boundaryCondition, forceIds);
///  y����������ڱ߽�
    lattice.periodicity().toggle(1,true);
    // Main loop over time iterations.
	plb_ofstream outfile("drag.dat");
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
        // At this point, the state of the lattice corresponds to the
        //   discrete time iT. However, the stored averages (getStoredAverageEnergy
        //   and getStoredAverageDensity) correspond to the previous time iT-1.

        if (iT%parameters.nStep(vtkSave)==0 && iT>0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; t=" << iT*parameters.getDeltaT();
        }
        
        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();
        T new_Energy = getStoredAverageEnergy(lattice);
        T Error_iT = std::fabs(old_Energy-new_Energy) / old_Energy;
        // At this point, the state of the lattice corresponds to the
        //   discrete time iT+1, and the stored averages are upgraded to time iT.
        if (iT%parameters.nStep(logT)==0) {
			      pcout << "; av energy ="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho ="
                  << getStoredAverageDensity<T>(lattice)
				  << endl;
				  T drag = lattice.getInternalStatistics().getSum(forceIds[0]);
                  T lift = lattice.getInternalStatistics().getSum(forceIds[1]);
                  T diameter = parameters.getResolution() + 1;
                  T vel = parameters.getLatticeU();
                  T avgRho = computeAverage(*computeDensity(lattice));
			      T drag_coef_factor = 1.0 / (0.5 * avgRho * util::sqr(vel) * diameter );
				  pcout << "Relative difference of Energy: " << setprecision(3)
                  << Error_iT <<std::endl;
                  pcout << "**********************************************" <<std::endl;
				  //����ļ���
				  outfile
				  << iT*parameters.getDeltaT()
				  << "; " << drag * drag_coef_factor
                  << "; " << lift * drag_coef_factor
				  << endl;
        }
		// ͨ��ƽ�����ܣ��ж��Ƿ�������
		if ( Error_iT < Tolerror )  {  
		pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);
		break;
		}
		old_Energy = new_Energy;
    }
    delete boundaryCondition;
}
