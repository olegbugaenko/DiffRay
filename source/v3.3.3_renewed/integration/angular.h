#pragma once
#include "matrix.h"
#include "../app.h"
#include "../lines.h"
#include <algorithm>
#include "../output/debugger.h"

struct AngleStep {
	double phi, theta, dphi, dtheta, Fi, dS;
	bool isPredicted;
	AngleStep(double cphi, double ctheta, double cdphi, double cdtheta)
	{
		phi = cphi;
		theta = ctheta;
		dphi = cdphi;
		dtheta = cdtheta;
		dS = cdphi * cdtheta;
		Fi = 0.;
		isPredicted = false;
	}

	AngleStep& operator=(const AngleStep& a)
	{
	    phi=a.phi;
	    theta=a.theta;
	    dphi=a.dphi;
	    dtheta=a.dtheta;
	    Fi = a.Fi;
	    isPredicted = a.isPredicted;
	    return *this;
	}
	
	double left()
	{
		return phi - dphi/2;
	}

	double right()
	{
		return phi + dphi/2;
	}

	double top()
	{
		return theta + dtheta/2;
	}

	double bottom()
	{
		return theta - dtheta/2;
	}

	bool contains(double phi, double theta)
	{
		if(phi < this->left()) return false;
		if(phi > this->right()) return false;
		if(theta > this->top()) return false;
		if(theta < this->bottom()) return false;
		return false;
	}
	double intersectionSquare(double phi2, double theta2, double dphi2, double dtheta2)
	{
		double left = std::max(this->left(), phi2 - dphi2/2);
		double bottom = std::max(this->bottom(), theta2 - dtheta2/2);
		double right = std::min(this->right(), phi2 + dphi2/2);
		double top = std::min(this->top(), theta2 + dtheta2/2);
		
		double width = right - left;
		double height = top - bottom;

		if(width < 0)
		    return 0;

		if(height < 0)
		    return 0;

		return width * height;
	} 
};

class Angular {
    public:
	double phi, theta, dphi, dtheta, nPhots, STot;
	AngleStep *angles[1000000];
	int AngleCount;
	double Vcenter, Scenter, SMcenter;
	
	//params are actual aperture constraints here
	Angular(double cphi, double ctheta, double cdphi, double cdtheta, int nSteps, bool doIterationOverSource)
	{
		CDebugger::debug("Init angular obj %le %le; [%le;%le]\n", cphi, ctheta, cdphi, cdtheta);
		phi = cphi;
		theta = ctheta;
		dphi = cdphi;
		dtheta = cdtheta;
		AngleCount = 0;
		long double phiW = dphi/nSteps;
		long double thetaW = dtheta/nSteps;
		
		for(int ip = 0; ip < nSteps; ip++)
		{
			for(int it = 0; it < nSteps; it++)
			{
				CDebugger::debug("inserting points: %d\n", AngleCount);
				long double aphi = cphi + (ip-(nSteps-1)/2)*phiW;
				long double atheta = ctheta + (it-(nSteps-1)/2)*thetaW;
				insertPoint(aphi, atheta, phiW, thetaW);
			}
		}
		if(App::usePredictiveMode)
		{
			//Run matrixes in predictive mode
			CDebugger::debug("Predictive\n");
			int nDepth = 0;
			while(predictiveIteration(doIterationOverSource) && nDepth < 3)
				nDepth++;
			sortByFi();
		}
	}

	bool sortByFi()
	{
		for(int i = 0; i < AngleCount; i++)
		{
			for(int j = 0; j < AngleCount-1; j++)
			{
				if(angles[j]->Fi > angles[j+1]->Fi)
				{
					AngleStep* buffer = angles[j];
					angles[j] = angles[j+1];
					angles[j+1] = buffer;
				}
			}
		}
		return true;
	}

	bool predictiveIteration(bool doIterationOverSource)
	{
		CDebugger::debug("start\n");
		for(int i = AngleCount - 1; i >= 0; i--)
		{
			if(angles[i]->isPredicted)
				continue;

			double currentPhi = angles[i]->phi;
			double currentTheta = angles[i]->theta;
			double currentDPhi = angles[i]->dphi;
			double currentDTheta = angles[i]->dtheta;
			CDebugger::debug(" -----===--- %d %le;%le\n",i,currentPhi, currentTheta);
			//now should setup matrix
			CMatrix::setup(currentPhi, currentTheta, currentDPhi, currentDTheta);
			CMatrix::run(false, doIterationOverSource);
			angles[i]->Fi = CMatrix::FFac/CMatrix::dS;
			angles[i]->isPredicted = true;
			CMatrix::freeMem();
		
		}
		sortByFi(); 
		double average = 0.;
		for(int ia=0;ia < AngleCount; ia++)
		{
			CDebugger::debug("Sorted: %le %le %le\n", angles[ia]->phi, angles[ia]->theta, angles[ia]->Fi);
			average += pow(angles[ia]->Fi,2.0);
		}
		//Get average
		average = sqrt(average/AngleCount);
		int itA = AngleCount - 1;
		while(itA > AngleCount/2 && (angles[itA]->Fi/average > AngleCount * 1./5))
		{
			//printf("Av: %le; Fi: %le\n",average,angles[itA]->Fi);
			itA--;
		}
		for(int ia = itA + 1; ia < AngleCount; ia++)
		{
			double currentPhi = angles[ia]->phi;
			double currentTheta = angles[ia]->theta;
			double currentDPhi = angles[ia]->dphi;
			double currentDTheta = angles[ia]->dtheta;
			removePoint(ia);
			int nSteps = 3;
			long double phiW = currentDPhi/nSteps;
			long double thetaW = currentDTheta/nSteps;
		
			for(int ip = 0; ip < nSteps; ip++)
			{
				for(int it = 0; it < nSteps; it++)
				{
					long double aphi = currentPhi + (ip-(nSteps-1)/2)*phiW;
					long double atheta = currentTheta + (it-(nSteps-1)/2)*thetaW;
					insertPoint(aphi, atheta, phiW, thetaW);
				}
			}
			CDebugger::log("Finished predictive iterations");
			return true;
		}
		return false;
	}

	void insertPoint(double cphi, double ctheta, double cdphi, double cdtheta)
	{
		angles[AngleCount] = new AngleStep(cphi,ctheta,cdphi,cdtheta);
		AngleCount++;
	}

	void removeLastPoint()
	{
		delete angles[AngleCount-1];
		AngleCount--;
	}
	
	void removePoint(int ip)
	{
		delete angles[ip];
		for(int i=ip;i<AngleCount-1;i++)
		    angles[i] = angles[i+1];
		
		angles[AngleCount-1] = NULL;
		AngleCount--;
	}

	void iterate(bool doIterationOverSource)
	{
		int acceptedCnt = 0;
		nPhots = 0.;
		STot = 0.;
		Vcenter = 0.;
		Scenter = 0.;
		SMcenter = 0.;
		CDebugger::log("Starting iterations");
		while(AngleCount > 0)
		{
			CDebugger::debug("AngleCount: %d (%d) ST: %d; nPhots: %le\n",AngleCount, acceptedCnt, CLine::iStat, nPhots);
			double currentPhi = angles[AngleCount-1]->phi;
			double currentTheta = angles[AngleCount-1]->theta;
			double currentDPhi = angles[AngleCount-1]->dphi;
			double currentDTheta = angles[AngleCount-1]->dtheta;
			removeLastPoint();
			//now should setup matrix
			CMatrix::setup(currentPhi, currentTheta, currentDPhi, currentDTheta);
			bool bWithStat = false;
			if(fabs(phi - currentPhi) < currentDPhi && fabs(theta - currentTheta) < currentDTheta)
				bWithStat = true;

			if(fabs(phi + dphi/2 - currentPhi) < currentDPhi &&
			    fabs(theta + dtheta/2 - currentTheta) < currentDTheta)
				bWithStat = true;
			if(fabs(phi - dphi/2 - currentPhi) < currentDPhi &&
			    fabs(theta - dtheta/2 - currentTheta) < currentDTheta)
				bWithStat = true;

			if(fabs(phi + dphi/2 - currentPhi) < currentDPhi &&
			    fabs(theta - dtheta/2 - currentTheta) < currentDTheta)
				bWithStat = true;
			if(fabs(phi - dphi/2 - currentPhi) < currentDPhi &&
			    fabs(theta + dtheta/2 - currentTheta) < currentDTheta)
				bWithStat = true;

			CMatrix::run(bWithStat, doIterationOverSource);
			bool bPassed = CMatrix::getDelta(nPhots);
			if(bPassed)
			{
				//add everything to total fluxes
				CMatrix::addFluxes(nPhots);
				STot += CMatrix::dS;
				acceptedCnt++;
				if(CMatrix::Vcenter > 1.e-50)
				{
					Vcenter += CMatrix::Vcenter;
					Scenter += CMatrix::Scenter;
					SMcenter += CMatrix::dS*App::distance*App::distance;
				}
			}
			else
			{
				for(int ip = 0; ip < 9; ip++)
				{
					insertPoint(CMatrix::rays[ip]->phi,
						CMatrix::rays[ip]->theta,
						CMatrix::rays[ip]->dphi,
						CMatrix::rays[ip]->dtheta
					);
				}
			}
			CMatrix::freeMem();
		}
	}
};