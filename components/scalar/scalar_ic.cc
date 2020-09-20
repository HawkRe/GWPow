#include "../../cosmo_includes.h"
#include "scalar_ic.h"
#include "SAMRAI/tbox/HDFDatabase.h"
#include "../elliptic_solver/full_multigrid.h"
#include "../elliptic_solver/multigrid_bd_handler.h"
#include "../../utils/Array.h"
#include <iostream>
#include <fstream>

#include "../horizon/AHFD/jtutil/util.hh"
#include "../horizon/AHFD/jtutil/array.hh"
#include "../horizon/AHFD/jtutil/util_String.h"
#include "../horizon/AHFD/jtutil/util_Table.h"
#include "../horizon/AHFD/jtutil/interpolator/InterpLocalUniform.h"
//#include "../horizon/AHFD/AHFD_types.h"
#include "../horizon/AHFD/AHFD_macros.h"

using namespace SAMRAI;

namespace cosmo
{

void scalar_ic_set_semianalytic_test(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_scalar_db)
{

}

  /*
    test for a spherical gaussian field collapse under 
    sommerfield boundary condition
  */
bool scalar_ic_set_scalar_collapse_sommerfield(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_scalar_db)
{
  return true;


}

bool scalar_ic_set_periodic_fast_collapse_test(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_scalar_db)
{

}

bool scalar_ic_set_scalar_collapse(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_scalar_db)
{
  return 0;
}


bool scalar_ic_set_oscillon(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_scalar_db)
{
      const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;


  const double * domain_lower = &grid_geometry.getXLower()[0];
  const double * domain_upper = &grid_geometry.getXUpper()[0];

  real_t L[3];
  double dx[3];
  
  for(int i = 0 ; i < DIM; i++)
  {
    L[i] = domain_upper[i] - domain_lower[i];
    dx[i] = (grid_geometry.getDx()[i]) / (1<<ln);
  }

  std::string boundary_type = "periodic";
  multigridBdHandler * bd_handler = new multigridBdHandler(boundary_type, L, 10);

  
  idx_t NX = round(L[0] / dx[0]);
  idx_t NY = round(L[1] / dx[1]); 
  idx_t NZ = round(L[2] / dx[2]);

  
  /******getting some parameters from input database****/

  real_t phi_0 = cosmo_scalar_db->getDoubleWithDefault("phi_0", 1.0);

  int n_max = cosmo_scalar_db->getIntegerWithDefault("n_max", 1);
  real_t delta_phi = 0.1;
  if(cosmo_scalar_db->keyExists("delta_phi"))
    delta_phi = cosmo_scalar_db->getDouble("delta_phi");
  real_t delta_phi_x = cosmo_scalar_db->getDoubleWithDefault("delta_phi_x", 0.1);
  real_t delta_phi_y = cosmo_scalar_db->getDoubleWithDefault("delta_phi_y", 0.1);
  real_t delta_phi_z = cosmo_scalar_db->getDoubleWithDefault("delta_phi_z", 0.1);


  /******ending collecting parameters from input database****/

  //  double * phi = new double[(NX+2*STENCIL_ORDER) * (NY+2*STENCIL_ORDER) * (NZ+2*STENCIL_ORDER)];
  CosmoArray<idx_t, real_t>  phi;
  phi.init(NX, NY, NZ);
  //  std::vector<double> phi(NX*NY*NZ, phi_0);
  
  LOOP3()
    phi[INDEX(i, j, k)] = phi_0;
/*
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<real_t> dist(0, 2.0*PI);
*/

  std::string initial_type =
    cosmo_scalar_db->getStringWithDefault("initial_type", "modes");
  if(initial_type == "modes")
  {
    for(int n = -n_max; n <= n_max; ++n)
    {
      if(n != 0)
      {
        // random phases
        /*
        real_t x_phase = dist(gen);
        real_t y_phase = dist(gen);
        real_t z_phase = dist(gen);
        std::ofstream random_phase;
	random_phase.open("random_phase.dat");
	random_phase << "This is added random phase: \n" <<  x_phase << '\n' << y_phase << '\n' << z_phase << '\n';
        random_phase.close();
	*/
	//given real number
	// real_t x_phase = 0.546719;
        // real_t y_phase = 0.726476;
        // real_t z_phase = 1.89999;
        real_t x_phase = PI;
        real_t y_phase = PI;
        real_t z_phase = PI;

        if(cosmo_scalar_db->keyExists("delta_phi"))
        {
          LOOP3()
          {
            // some sinusoidal modes
            phi[INDEX(i,j,k)] += delta_phi*(
              cos(2.0*PI*((real_t) n/NX)*((real_t)i + 0.5) + x_phase)
              + cos(2.0*PI*((real_t) n/NY)*((real_t)j + 0.5) + y_phase)
              + cos(2.0*PI*((real_t) n/NZ)*((real_t)k + 0.5) + z_phase)
            );
          }
        }
        else
        {
          LOOP3()
          {
          phi[INDEX(i,j,k)] +=
            delta_phi_x * cos(2.0*PI*((real_t) n/NX)*((real_t)i + 0.5) + x_phase)
            + delta_phi_y * cos(2.0*PI*((real_t) n/NY)*((real_t)j + 0.5) + y_phase)
            + delta_phi_z * cos(2.0*PI*((real_t) n/NZ)*((real_t)k + 0.5) + z_phase);
          }
        }
      }
    }
  }
  else if(initial_type == "gaussian")
  {
    double r0 = cosmo_scalar_db->getDoubleWithDefault("r0", 0);
    double sigma = cosmo_scalar_db->getDoubleWithDefault("sigma", 1.0);
    double q = cosmo_scalar_db->getDoubleWithDefault("q", 2.0);
    LOOP3()
    {
      double x = L[0] / NX * ((double)i + 0.5) - L[0] / 2.0;
      double y = L[1] / NX * ((double)j + 0.5) - L[1] / 2.0;
      double z = L[2] / NX * ((double)k + 0.5) - L[2] / 2.0;
      double r = sqrt(x * x + y * y + z * z);
      phi[INDEX(i,j,k)] = delta_phi * pw3(r) *
      exp( - pow(fabs( (r - r0) / sigma) , q)) ;
    }
  }
    


  bd_handler->fillBoundary(phi._array, phi.nx, phi.ny, phi.nz);

  double tot_r = 0, tot_v = 0.0;
  real_t rho_sigma = 0;
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
    {
      const std::shared_ptr<hier::Patch> & patch = *pit;

      bssn->initPData(patch);
      bssn->initMDA(patch);

      scalar->initPData(patch);
      scalar->initMDA(patch);
    
      arr_t & phi_a = scalar->phi_a;
      arr_t & a_a = bssn->a_a; // field
    
      const hier::Box& box = bssn->DIFFchi_a_pdata->getGhostBox();
      const hier::Box& inner_box = patch->getBox();

      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];

      const int * inner_lower = &inner_box.lower()[0];
      const int * inner_upper = &inner_box.upper()[0];

    
      for(int k = lower[2]; k <= upper[2]; k++)
	{
	  for(int j = lower[1]; j <= upper[1]; j++)
	    {
	      for(int i = lower[0]; i <= upper[0]; i++)
		{

                  phi_a(i, j, k) = phi[INDEX(i,j,k)];
                  a_a(i, j, k) = 1.0;
		}
	    }
	}

      for(int k = inner_lower[2]; k <= inner_upper[2]; k++)
	{
	  for(int j = inner_lower[1]; j <= inner_upper[1]; j++)
	    {
	      for(int i = inner_lower[0]; i <= inner_upper[0]; i++)
		{         
		  BSSNData bd = {0};
		  ScalarData sd = {0};

		  sd.phi = phi[INDEX(i,j,k)];
		  tot_r += 
		    (0.5 *  (
							  (pw2((1.0/12.0*phi[INDEX(i-2,j,k)] - 2.0/3.0*phi[INDEX(i-1,j,k)] + 2.0/3.0*phi[INDEX(i+1,j,k)]- 1.0/12.0*phi[INDEX(i+2,j,k)])/dx[0])
							   + pw2((1.0/12.0*phi[INDEX(i,j-2,k)] - 2.0/3.0*phi[INDEX(i,j-1,k)] + 2.0/3.0*phi[INDEX(i,j+1,k)]- 1.0/12.0*phi[INDEX(i,j+2,k)])/dx[1])
							   +pw2((1.0/12.0*phi[INDEX(i,j,k-2)] - 2.0/3.0*phi[INDEX(i,j,k-1)] + 2.0/3.0*phi[INDEX(i,j,k+1)]- 1.0/12.0*phi[INDEX(i,j,k+2)])/dx[2]))
							  )
		     + scalar->potentialHandler->ev_potential(&bd, &sd));

		}
	    }
	}   
    }

  mpi.AllReduce(&tot_r,1,MPI_SUM);

  double avg_r = tot_r / ((double)NX * NY * NZ);
  

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
    {
      const std::shared_ptr<hier::Patch> & patch = *pit;
 
      bssn->initPData(patch);
      bssn->initMDA(patch);
 
      scalar->initPData(patch);
      scalar->initMDA(patch);
     
      arr_t & phi_a = scalar->phi_a; // field
      arr_t & H_a = bssn->H_a; // field``
      
      const hier::Box& box = bssn->DIFFchi_a_pdata->getGhostBox();
      const hier::Box& inner_box = patch->getBox();
 
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];
 
      const int * inner_lower = &inner_box.lower()[0];
      const int * inner_upper = &inner_box.upper()[0];
 
     
      for(int k = inner_lower[2]; k <= inner_upper[2]; k++)
	{
	  for(int j = inner_lower[1]; j <= inner_upper[1]; j++)
	    {
	      for(int i = inner_lower[0]; i <= inner_upper[0]; i++)
		{
		  BSSNData bd = {0};
		  ScalarData sd = {0};
                  H_a(i, j, k) = sqrt(8.0 * PI * avg_r/3.0);
		  sd.phi = phi[INDEX(i,j,k)];
 
		  rho_sigma += 
		    pw2((0.5  * (
							      (pw2((1.0/12.0*phi[INDEX(i-2,j,k)] - 2.0/3.0*phi[INDEX(i-1,j,k)] + 2.0/3.0*phi[INDEX(i+1,j,k)]- 1.0/12.0*phi[INDEX(i+2,j,k)])/dx[0])
							       + pw2((1.0/12.0*phi[INDEX(i,j-2,k)] - 2.0/3.0*phi[INDEX(i,j-1,k)] + 2.0/3.0*phi[INDEX(i,j+1,k)]- 1.0/12.0*phi[INDEX(i,j+2,k)])/dx[1])
							       +pw2((1.0/12.0*phi[INDEX(i,j,k-2)] - 2.0/3.0*phi[INDEX(i,j,k-1)] + 2.0/3.0*phi[INDEX(i,j,k+1)]- 1.0/12.0*phi[INDEX(i,j,k+2)])/dx[2]))
							      )
			 + scalar->potentialHandler->ev_potential(&bd, &sd)) - avg_r);
		}
	    }
	}
    }
 
  mpi.AllReduce(&rho_sigma,1,MPI_SUM);
  tbox::pout<<"sigma rho is "
            <<sqrt(pw3(NX) / (pw3(NX)-1) * rho_sigma / ((double)NX * NY * NZ))<<
    " sigma rho / rho is "<<sqrt(pw3(NX) / (pw3(NX)-1) * rho_sigma / ((double)NX * NY * NZ)) /
    ( avg_r)<<"\n";

  
  return true;



























  
  // double tot_r = 0, tot_v = 0.0;
  // real_t rho_sigma = 0.0;
  
  // for( hier::PatchLevel::iterator pit(level->begin());
  //      pit != level->end(); ++pit)
  // {
  //   const std::shared_ptr<hier::Patch> & patch = *pit;

  //   bssn->initPData(patch);
  //   bssn->initMDA(patch);

  //   scalar->initPData(patch);
  //   scalar->initMDA(patch);

  //   arr_t & phi_a = scalar->phi_a; // field
  //   arr_t & a_a = bssn->a_a; // field
    
  //   const hier::Box& box = bssn->DIFFchi_a_pdata->getGhostBox();
  //   const hier::Box& inner_box = patch->getBox();

  //   const int * lower = &box.lower()[0];
  //   const int * upper = &box.upper()[0];

  //   const int * inner_lower = &inner_box.lower()[0];
  //   const int * inner_upper = &inner_box.upper()[0];

  //   for(int k = lower[2]; k <= upper[2]; k++)
  //   {
  //     for(int j = lower[1]; j <= upper[1]; j++)
  //     {
  //       for(int i = lower[0]; i <= upper[0]; i++)
  //       {
  //         phi_a(i, j, k) = phi[INDEX(i,j,k)];
  //         a_a(i, j, k) = 1.0;
  //         // if(i == 129 && j == 0 && k == 0)
  //         //   std::cout<<"phi_a(129, 0, 0) "<<phi_a(129, 0, 0)<<" "<<r_field[fft_index]<<" "<<fft_index<<"\n";
  //       }
  //     }
  //   }
  //   //    std::cout<<"phi_a(129, 0, 0) "<<phi_a(129, 0, 0)<<"\n";
  //   for(int k = inner_lower[2]; k <= inner_upper[2]; k++)
  //   {
  //     for(int j = inner_lower[1]; j <= inner_upper[1]; j++)
  //     {
  //       for(int i = inner_lower[0]; i <= inner_upper[0]; i++)
  //       {
  //         BSSNData bd = {0};
  //         ScalarData sd = {0};

  //         sd.phi = phi_a(i,j,k);
  //         tot_r += 0.5*
  //           (
  //             (pw2((1.0/12.0*phi_a(i-2,j,k) - 2.0/3.0*phi_a(i-1,j,k) + 2.0/3.0*phi_a(i+1,j,k)- 1.0/12.0*phi_a(i+2,j,k))/dx[0])
  //              + pw2((1.0/12.0*phi_a(i,j-2,k) - 2.0/3.0*phi_a(i,j-1,k) + 2.0/3.0*phi_a(i,j+1,k)- 1.0/12.0*phi_a(i,j+2,k))/dx[1])
  //              +pw2((1.0/12.0*phi_a(i,j,k-2) - 2.0/3.0*phi_a(i,j,k-1) + 2.0/3.0*phi_a(i,j,k+1)- 1.0/12.0*phi_a(i,j,k+2))/dx[2]))
  //           ) + scalar->potentialHandler->ev_potential(&bd, &sd);

  //         // if(k == 0 && j == 0 )
  //         //   std::cout<<" "<<box<<" "<<inner_box<<" "<<i<<" "<<j<<" "<<k<<" "<<tot_r<<" "<<            (
  //         //     (pw2(- 1.0/12.0*phi_a(i+2,j,k))
  //         //      )
  //         //   )
  //         //          <<" "<<scalar->potentialHandler->ev_potential(&bd, &sd)<<"  "<<sd.phi<<"\n";
  //       }
  //     }
  //   }    
  // }

  // mpi.AllReduce(&tot_r,1,MPI_SUM);
  // tot_r = tot_r / ((double)NX * NY * NZ);
  // std::cout<<tot_r<<"\n";

  // for( hier::PatchLevel::iterator pit(level->begin());
  //      pit != level->end(); ++pit)
  // {
  //   const std::shared_ptr<hier::Patch> & patch = *pit;

  //   bssn->initPData(patch);
  //   bssn->initMDA(patch);

  //   scalar->initPData(patch);
  //   scalar->initMDA(patch);

  //   arr_t & H_a = bssn->H_a; // field
    
  //   const hier::Box& box = bssn->DIFFchi_a_pdata->getGhostBox();
  //   const hier::Box& inner_box = patch->getBox();

  //   const int * lower = &box.lower()[0];
  //   const int * upper = &box.upper()[0];

  //   const int * inner_lower = &inner_box.lower()[0];
  //   const int * inner_upper = &inner_box.upper()[0];

  //   for(int k = lower[2]; k <= upper[2]; k++)
  //   {
  //     for(int j = lower[1]; j <= upper[1]; j++)
  //     {
  //       for(int i = lower[0]; i <= upper[0]; i++)
  //       {
  //         H_a(i, j, k) = sqrt(8.0 * PI * tot_r/3.0);
  //       }
  //     }
  //   }
  // }  

  // return true;
}

bool scalar_ic_set_scalar_gaussian_random(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_scalar_db)
{
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;


  const double * domain_lower = &grid_geometry.getXLower()[0];
  const double * domain_upper = &grid_geometry.getXUpper()[0];

  real_t L[3];
  double dx[3];
  
  for(int i = 0 ; i < DIM; i++)
  {
    L[i] = domain_upper[i] - domain_lower[i];
    dx[i] = (grid_geometry.getDx()[i]) / (1<<ln);
  }

  std::string boundary_type = "periodic";
  multigridBdHandler * bd_handler = new multigridBdHandler(boundary_type, L, 10);
  
  idx_t NX = round(L[0] / dx[0]);
  idx_t NY = round(L[1] / dx[1]); 
  idx_t NZ = round(L[2] / dx[2]);

  /******getting some parameters from input database****/

  real_t rho_K_matter = 3.0/PI/8.0;
  real_t peak_amplitude_frac = cosmo_scalar_db->getDoubleWithDefault("peak_amplitude_frac", 0.0);
  real_t peak_amplitude = peak_amplitude_frac*(1.0); // scaling in arb. units
  real_t q_coef = cosmo_scalar_db->getDoubleWithDefault("q_coef", 0.0);
  real_t ic_spec_cut = cosmo_scalar_db->getDoubleWithDefault("ic_spec_cut", 0.0);

  real_t peak_k = 100;

  real_t phi_0 = cosmo_scalar_db->getDoubleWithDefault("phi_0", 1.0);

  int num_vcycles = cosmo_scalar_db->getIntegerWithDefault("vcycles", 20);

  real_t relaxation_tolerance = cosmo_scalar_db->getDoubleWithDefault("relaxation_tolerance", 1e-8);

  double DIFFalpha_0 = cosmo_scalar_db->getDoubleWithDefault("DIFFalpha", 0);
  
  /*********Finishing getting input parameters*******/

  real_t px, py, pz, pmag;
  real_t scale;

  
  // populate "field" with random values
  std::random_device rd;
  std::mt19937 gen(9.0 /* rd()*/);
  std::normal_distribution<real_t> gaussian_distribution;
  std::uniform_real_distribution<double> angular_distribution(0.0, 2.0*PI);
  // calling these here before looping suppresses a warning (bug)
  gaussian_distribution(gen);
  angular_distribution(gen);

  fftw_complex *f_field;

  double *r_field = new double[NX*NY*NZ];

  memset(r_field, 0, NX*NY*NZ*sizeof(double));

  
  f_field = (fftw_complex *) fftw_malloc(NX*NY*(NZ/2+1)
                                         *((long long) sizeof(fftw_complex)));
  // plans for taking FFTs
  fftw_plan p_c2r;

  p_c2r = fftw_plan_dft_c2r_3d(NX, NY, NZ,
                               f_field, r_field,
                               FFTW_MEASURE);

  // scale amplitudes in fourier space
  // don't expect to run at >512^3 anytime soon; loop over all momenta out to that resolution.
  // this won't work for a larger grid.
  idx_t NMAX = 512;
  real_t max_fft_index = NX*NY*(NZ/2+1) - 1;
  for(int i=0; i<NMAX; i++)
  {
    px = (real_t) (i<=NMAX/2 ? i : i-NMAX);
    for(int j=0; j<NMAX; j++)
    {
      py = (real_t) (j<=NMAX/2 ? j : j-NMAX);
      for(int k=0; k<NMAX/2+1; k++)
      {
        pz = (real_t) k;

        // generate the same random modes for all resolutions (up to NMAX)
        //        real_t rand_mag = gaussian_distribution(gen);
        real_t rand_phase = angular_distribution(gen);

        // only store momentum values for relevant bins
        if( fabs(px) < (real_t) NX/2+1 + 0.01
            && fabs(py) < (real_t) NY/2+1 + 0.01
            && fabs(pz) < (real_t) NZ/2+1 + 0.01)
        {
          idx_t fft_index = FFT_NP_INDEX(
            px > -0.5 ? ROUND_2_IDXT(px) : (NX + ROUND_2_IDXT(px)),
            py > -0.5 ? ROUND_2_IDXT(py) : (NY + ROUND_2_IDXT(py)),
            pz > -0.5 ? ROUND_2_IDXT(pz) : (NZ + ROUND_2_IDXT(pz))
          );
 
          if(fft_index > max_fft_index)
          {
            tbox::pout<< "Warning: index " << fft_index << " is greater than max ("
                      << max_fft_index << ").\n";
            // std::cout<<px<<" "<<py<<" "<<pz<<"\n";
            // std::cout<<(px > -0.5 ? ROUND_2_IDXT(px) : (NX + ROUND_2_IDXT(px)))
            //          <<" "<<(py > -0.5 ? ROUND_2_IDXT(py) : (NY + ROUND_2_IDXT(py)))
            //          <<" "<<(pz > -0.5 ? ROUND_2_IDXT(pz) : (NZ + ROUND_2_IDXT(pz)))<<"\n";
            fft_index = max_fft_index;
          }

          pmag = sqrt(
            pw2(px * ( (real_t) NX / (real_t) NX ) )
            + pw2(py * ( (real_t) NX / (real_t) NY ))
            + pw2(pz * ( (real_t) NX / (real_t) NZ ))
          );
          // Scale by power spectrum
          // don't want much power on scales smaller than ~3 pixels
          // Or scales p > 1/(3*dx)
          real_t cutoff = 1.0 / (
            1.0 + exp(10.0*(pmag - ic_spec_cut))
          );
          real_t pre = peak_amplitude;
          //real_t cosmo_power_spectrum =  pre/pow(fabs(pmag)/peak_k, 3.0);
          //real_t cosmo_power_spectrum =  pre/(1.0 + pow(fabs(pmag)/peak_k, 4.0)/3.0)/pow(fabs(pmag)/peak_k, 3.0);
          //real_t cosmo_power_spectrum =  pre;
          real_t cosmo_power_spectrum =  pw3(1) / (2.0 * sqrt(q_coef + pw2((2.0 * PI / L[0]) * pmag))) ;
          //          real_t cosmo_power_spectrum = 1.0;
          scale = cutoff;
          
          std::weibull_distribution<> d(2, sqrt(2.0 *cosmo_power_spectrum) );
          //real_t rand_mag = d(gen);
          real_t rand_mag = gaussian_distribution(gen) * sqrt(cosmo_power_spectrum);
          // if(pmag <= 2)
          //   rand_mag = sqrt(cosmo_power_spectrum) *1.5;
          

          f_field[fft_index][0] = scale*rand_mag*cos(rand_phase);
          f_field[fft_index][1] = scale*rand_mag*sin(rand_phase);


         
        }
	else if( fabs(px) < (real_t) 256/2+1 + 0.01 && fabs(py) < (real_t) 256/2+1 + 0.01 && fabs(pz) < (real_t) 256/2+1 + 0.01)
	  real_t rand_mag = gaussian_distribution(gen);
      }
    }
  }

  for(int i = 0; i < NX; i++)
    std::cout<<(f_field)[FFT_NP_INDEX(i,60,60)][0]<<" ";
  std::cout<<"\n";
  // zero-mode (mean density)... set this to something later
  (f_field)[FFT_NP_INDEX(0,0,0)][0] = 0;
  (f_field)[FFT_NP_INDEX(0,0,0)][1] = 0;

  
  fftw_execute_dft_c2r(p_c2r, f_field, r_field);

  // for(int i = 0; i < NX;i ++)
  // {
  //   idx_t fft_index = NP_INDEX(i,0,0);
  //   std::cout<<r_field[fft_index]<<" ";
  // }
  // std::cout<<"\n";

  CosmoArray<idx_t, real_t>  phi;
  phi.init(NX, NY, NZ);
  double max_phi = -1;
  LOOP3()
  {
    idx_t fft_index = NP_INDEX(i,j,k);
    phi[INDEX(i, j, k)] = phi_0 + 16.0*r_field[fft_index]/pow(L[0],1.5);
    // if(i == 0 && j == 0)
    //   std::cout<<r_field[fft_index]<<" ";
    // if(i == 0 && j == 0)
    //   std::cout<<"\n";

    max_phi = std::max(max_phi, phi[INDEX(i, j, k)]);
  }
  std::cout<<phi_0<<" "<<max_phi<<"\n";
  bd_handler->fillBoundary(phi._array, phi.nx, phi.ny, phi.nz);
  
  double tot_r = 0, tot_v = 0.0;
  real_t rho_sigma = 0.0;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    bssn->initPData(patch);
    bssn->initMDA(patch);

    scalar->initPData(patch);
    scalar->initMDA(patch);

    arr_t & phi_a = scalar->phi_a; // field
    arr_t & a_a = bssn->a_a; // field
    
    const hier::Box& box = bssn->DIFFchi_a_pdata->getGhostBox();
    const hier::Box& inner_box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    const int * inner_lower = &inner_box.lower()[0];
    const int * inner_upper = &inner_box.upper()[0];

    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          phi_a(i, j, k) = phi[INDEX(i,j,k)];
          a_a(i, j, k) = 1.0;
          // if(i == 129 && j == 0 && k == 0)
          //   std::cout<<"phi_a(129, 0, 0) "<<phi_a(129, 0, 0)<<" "<<r_field[fft_index]<<" "<<fft_index<<"\n";
        }
      }
    }
    //    std::cout<<"phi_a(129, 0, 0) "<<phi_a(129, 0, 0)<<"\n";
    for(int k = inner_lower[2]; k <= inner_upper[2]; k++)
    {
      for(int j = inner_lower[1]; j <= inner_upper[1]; j++)
      {
        for(int i = inner_lower[0]; i <= inner_upper[0]; i++)
        {
          BSSNData bd = {0};
          ScalarData sd = {0};

          sd.phi = phi_a(i,j,k);
          tot_r += 0.5*
            (
              (pw2((1.0/12.0*phi_a(i-2,j,k) - 2.0/3.0*phi_a(i-1,j,k) + 2.0/3.0*phi_a(i+1,j,k)- 1.0/12.0*phi_a(i+2,j,k))/dx[0])
               + pw2((1.0/12.0*phi_a(i,j-2,k) - 2.0/3.0*phi_a(i,j-1,k) + 2.0/3.0*phi_a(i,j+1,k)- 1.0/12.0*phi_a(i,j+2,k))/dx[1])
               +pw2((1.0/12.0*phi_a(i,j,k-2) - 2.0/3.0*phi_a(i,j,k-1) + 2.0/3.0*phi_a(i,j,k+1)- 1.0/12.0*phi_a(i,j,k+2))/dx[2]))
            ) + scalar->potentialHandler->ev_potential(&bd, &sd);

          // if(k == 0 && j == 0 )
          //   std::cout<<" "<<box<<" "<<inner_box<<" "<<i<<" "<<j<<" "<<k<<" "<<tot_r<<" "<<            (
          //     (pw2(- 1.0/12.0*phi_a(i+2,j,k))
          //      )
          //   )
          //          <<" "<<scalar->potentialHandler->ev_potential(&bd, &sd)<<"  "<<sd.phi<<"\n";
        }
      }
    }    
  }

  mpi.AllReduce(&tot_r,1,MPI_SUM);
  tot_r = tot_r / ((double)NX * NY * NZ);
  std::cout<<tot_r<<"\n";

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    bssn->initPData(patch);
    bssn->initMDA(patch);

    scalar->initPData(patch);
    scalar->initMDA(patch);

    arr_t & H_a = bssn->H_a; // field
    
    const hier::Box& box = bssn->DIFFchi_a_pdata->getGhostBox();
    const hier::Box& inner_box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    const int * inner_lower = &inner_box.lower()[0];
    const int * inner_upper = &inner_box.upper()[0];

    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          H_a(i, j, k) = sqrt(8.0 * PI * tot_r/3.0);
        }
      }
    }
  }
  return 1;
}

bool scalar_ic_set_scalar_gaussian_random_defrost(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_scalar_db)
{
  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

  std::shared_ptr<hier::PatchLevel> level(
    hierarchy->getPatchLevel(ln));

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;


  const double * domain_lower = &grid_geometry.getXLower()[0];
  const double * domain_upper = &grid_geometry.getXUpper()[0];

  real_t L[3];
  double dx[3];
  
  for(int i = 0 ; i < DIM; i++)
  {
    L[i] = domain_upper[i] - domain_lower[i];
    dx[i] = (grid_geometry.getDx()[i]) / (1<<ln);
  }

  std::string boundary_type = "periodic";
  multigridBdHandler * bd_handler = new multigridBdHandler(boundary_type, L, 10);
  
  idx_t NX = round(L[0] / dx[0]);
  idx_t NY = round(L[1] / dx[1]); 
  idx_t NZ = round(L[2] / dx[2]);

  /******getting some parameters from input database****/

  real_t rho_K_matter = 3.0/PI/8.0;
  real_t peak_amplitude_frac = cosmo_scalar_db->getDoubleWithDefault("peak_amplitude_frac", 0.0);
  real_t peak_amplitude = peak_amplitude_frac*(1.0); // scaling in arb. units
  real_t q_coef = cosmo_scalar_db->getDoubleWithDefault("q_coef", 0.0);
  real_t ic_spec_cut = cosmo_scalar_db->getDoubleWithDefault("ic_spec_cut", 0.0);

  real_t peak_k = 100;

  real_t phi_0 = cosmo_scalar_db->getDoubleWithDefault("phi_0", 1.0);

  int num_vcycles = cosmo_scalar_db->getIntegerWithDefault("vcycles", 20);

  real_t relaxation_tolerance = cosmo_scalar_db->getDoubleWithDefault("relaxation_tolerance", 1e-8);

  double DIFFalpha_0 = cosmo_scalar_db->getDoubleWithDefault("DIFFalpha", 0);
  
  /*********Finishing getting input parameters*******/

  real_t px, py, pz, pmag;
  real_t scale;

  
  // populate "field" with random values
  std::random_device rd;
  std::mt19937 gen(9.0 /* rd()*/);
  std::normal_distribution<real_t> gaussian_distribution;
  std::uniform_real_distribution<double> angular_distribution(0.0, 2.0*PI);
  // calling these here before looping suppresses a warning (bug)
  gaussian_distribution(gen);
  angular_distribution(gen);




  /*****************begin defrost part***********************************/

  field_size fs(NX);
  field<double> phi_df;
  fs.calculate_size_totals();
  phi_df.construct(fs, false);

  initializer<double> *init =  (initializer<double> *) new defrost_style_initializer<double>(fs,
                                                          phi_df);

  BSSNData bd = {0};
  ScalarData sd = {0};

  sd.phi = phi_0;
  

  init->initialize(sqrt(q_coef), sqrt(8.0 * PI * scalar->potentialHandler->ev_potential(&bd, &sd)/3.0), L[0]);

  phi_df.mdata[0][0] =  fs.total_gridpoints * phi_0;
  phi_df.mdata[0][1] = 0.;

  
  phi_df.switch_state(position, true);

  
  for(int i = 0; i < 100; i++)
    std::cout<<phi_df.mdata[i][0]<<" ";
  std::cout<<"\n";
  
  std::cout<<"!!!!!!!!! "<<fs.total_gridpoints<<" "<<sqrt(8.0 * PI * scalar->potentialHandler->ev_potential(&bd, &sd)/3.0)<<" "
           <<phi_df.ldl<<"\n";
  
  /*****************end defronst part************************************/


  
  fftw_complex *f_field;

  double *r_field = new double[NX*NY*NZ];

  memset(r_field, 0, NX*NY*NZ*sizeof(double));

  
  f_field = (fftw_complex *) fftw_malloc(NX*NY*(NZ/2+1)
                                         *((long long) sizeof(fftw_complex)));
  // plans for taking FFTs
  fftw_plan p_c2r;

  p_c2r = fftw_plan_dft_c2r_3d(NX, NY, NZ,
                               f_field, r_field,
                               FFTW_MEASURE);


  
  
  // scale amplitudes in fourier space
  // don't expect to run at >512^3 anytime soon; loop over all momenta out to that resolution.
  // this won't work for a larger grid.
  idx_t NMAX = 512;
  real_t max_fft_index = NX*NY*(NZ/2+1) - 1;
  for(int i=0; i<NMAX; i++)
  {
    px = (real_t) (i<=NMAX/2 ? i : i-NMAX);
    for(int j=0; j<NMAX; j++)
    {
      py = (real_t) (j<=NMAX/2 ? j : j-NMAX);
      for(int k=0; k<NMAX/2+1; k++)
      {
        pz = (real_t) k;

        // generate the same random modes for all resolutions (up to NMAX)
        //        real_t rand_mag = gaussian_distribution(gen);
        real_t rand_phase = angular_distribution(gen);

        // only store momentum values for relevant bins
        if( fabs(px) < (real_t) NX/2+1 + 0.01
            && fabs(py) < (real_t) NY/2+1 + 0.01
            && fabs(pz) < (real_t) NZ/2+1 + 0.01)
        {
          idx_t fft_index = FFT_NP_INDEX(
            px > -0.5 ? ROUND_2_IDXT(px) : (NX + ROUND_2_IDXT(px)),
            py > -0.5 ? ROUND_2_IDXT(py) : (NY + ROUND_2_IDXT(py)),
            pz > -0.5 ? ROUND_2_IDXT(pz) : (NZ + ROUND_2_IDXT(pz))
          );
 
          if(fft_index > max_fft_index)
          {
            tbox::pout<< "Warning: index " << fft_index << " is greater than max ("
                      << max_fft_index << ").\n";
            // std::cout<<px<<" "<<py<<" "<<pz<<"\n";
            // std::cout<<(px > -0.5 ? ROUND_2_IDXT(px) : (NX + ROUND_2_IDXT(px)))
            //          <<" "<<(py > -0.5 ? ROUND_2_IDXT(py) : (NY + ROUND_2_IDXT(py)))
            //          <<" "<<(pz > -0.5 ? ROUND_2_IDXT(pz) : (NZ + ROUND_2_IDXT(pz)))<<"\n";
            fft_index = max_fft_index;
          }

          pmag = sqrt(
            pw2(px * ( (real_t) NX / (real_t) NX ) )
            + pw2(py * ( (real_t) NX / (real_t) NY ))
            + pw2(pz * ( (real_t) NX / (real_t) NZ ))
          );
          // Scale by power spectrum
          // don't want much power on scales smaller than ~3 pixels
          // Or scales p > 1/(3*dx)
          real_t cutoff = 1.0 / (
            1.0 + exp(10.0*(pmag - ic_spec_cut))
          );
          real_t pre = peak_amplitude;
          //real_t cosmo_power_spectrum =  pre/pow(fabs(pmag)/peak_k, 3.0);
          //real_t cosmo_power_spectrum =  pre/(1.0 + pow(fabs(pmag)/peak_k, 4.0)/3.0)/pow(fabs(pmag)/peak_k, 3.0);
          //real_t cosmo_power_spectrum =  pre;
          real_t cosmo_power_spectrum =  pw3(1) / (2.0 * sqrt(q_coef + pw2((2.0 * PI / L[0]) * pmag))) ;
          //          real_t cosmo_power_spectrum = 1.0;
          scale = cutoff;
          
          std::weibull_distribution<> d(2, sqrt(2.0 *cosmo_power_spectrum) );
          //real_t rand_mag = d(gen);
          real_t rand_mag = gaussian_distribution(gen) * sqrt(cosmo_power_spectrum);
          // if(pmag <= 2)
          //   rand_mag = sqrt(cosmo_power_spectrum) *1.5;
          

          f_field[fft_index][0] = scale*rand_mag*cos(rand_phase);
          f_field[fft_index][1] = scale*rand_mag*sin(rand_phase);


         
        }
      }
    }
  }

  for(int i = 0; i < NX; i++)
    std::cout<<(f_field)[FFT_NP_INDEX(i,60,60)][0]<<" ";
  std::cout<<"\n";
  // zero-mode (mean density)... set this to something later
  (f_field)[FFT_NP_INDEX(0,0,0)][0] = 0;
  (f_field)[FFT_NP_INDEX(0,0,0)][1] = 0;

  
  fftw_execute_dft_c2r(p_c2r, f_field, r_field);

  // for(int i = 0; i < NX;i ++)
  // {
  //   idx_t fft_index = NP_INDEX(i,0,0);
  //   std::cout<<r_field[fft_index]<<" ";
  // }
  // std::cout<<"\n";

  CosmoArray<idx_t, real_t>  phi;
  phi.init(NX, NY, NZ);
  double max_phi = -1;
  LOOP3()
  {
    idx_t fft_index = NP_INDEX(i,j,k);
    int fdx = k + phi_df.ldl*(j + fs.n*i);
    //       phi[INDEX(i, j, k)] = phi_0 + r_field[fft_index]/pw3(NX);
        phi[INDEX(i, j, k)] = phi_df.data[fdx];
    // if(i == 0 && j == 0)
    //   std::cout<<r_field[fft_index]<<" ";
    // if(i == 0 && j == 0)
    //   std::cout<<"\n";

    max_phi = std::max(max_phi, phi[INDEX(i, j, k)]);
  }
  std::cout<<phi_0<<" "<<max_phi<<"\n";
  bd_handler->fillBoundary(phi._array, phi.nx, phi.ny, phi.nz);
  
  double tot_r = 0, tot_v = 0.0;
  real_t rho_sigma = 0.0;
  
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    bssn->initPData(patch);
    bssn->initMDA(patch);

    scalar->initPData(patch);
    scalar->initMDA(patch);

    arr_t & phi_a = scalar->phi_a; // field
    arr_t & a_a = bssn->a_a; // field
    
    const hier::Box& box = bssn->DIFFchi_a_pdata->getGhostBox();
    const hier::Box& inner_box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    const int * inner_lower = &inner_box.lower()[0];
    const int * inner_upper = &inner_box.upper()[0];

    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          phi_a(i, j, k) = phi[INDEX(i,j,k)];
          a_a(i, j, k) = 1.0;
          // if(i == 129 && j == 0 && k == 0)
          //   std::cout<<"phi_a(129, 0, 0) "<<phi_a(129, 0, 0)<<" "<<r_field[fft_index]<<" "<<fft_index<<"\n";
        }
      }
    }
    //    std::cout<<"phi_a(129, 0, 0) "<<phi_a(129, 0, 0)<<"\n";
    for(int k = inner_lower[2]; k <= inner_upper[2]; k++)
    {
      for(int j = inner_lower[1]; j <= inner_upper[1]; j++)
      {
        for(int i = inner_lower[0]; i <= inner_upper[0]; i++)
        {
          BSSNData bd = {0};
          ScalarData sd = {0};

          sd.phi = phi_a(i,j,k);
          tot_r += 0.5*
            (
              (pw2((1.0/12.0*phi_a(i-2,j,k) - 2.0/3.0*phi_a(i-1,j,k) + 2.0/3.0*phi_a(i+1,j,k)- 1.0/12.0*phi_a(i+2,j,k))/dx[0])
               + pw2((1.0/12.0*phi_a(i,j-2,k) - 2.0/3.0*phi_a(i,j-1,k) + 2.0/3.0*phi_a(i,j+1,k)- 1.0/12.0*phi_a(i,j+2,k))/dx[1])
               +pw2((1.0/12.0*phi_a(i,j,k-2) - 2.0/3.0*phi_a(i,j,k-1) + 2.0/3.0*phi_a(i,j,k+1)- 1.0/12.0*phi_a(i,j,k+2))/dx[2]))
            ) + scalar->potentialHandler->ev_potential(&bd, &sd);

          // if(k == 0 && j == 0 )
          //   std::cout<<" "<<box<<" "<<inner_box<<" "<<i<<" "<<j<<" "<<k<<" "<<tot_r<<" "<<            (
          //     (pw2(- 1.0/12.0*phi_a(i+2,j,k))
          //      )
          //   )
          //          <<" "<<scalar->potentialHandler->ev_potential(&bd, &sd)<<"  "<<sd.phi<<"\n";
        }
      }
    }    
  }

  mpi.AllReduce(&tot_r,1,MPI_SUM);
  tot_r = tot_r / ((double)NX * NY * NZ);
  std::cout<<tot_r<<"\n";

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;

    bssn->initPData(patch);
    bssn->initMDA(patch);

    scalar->initPData(patch);
    scalar->initMDA(patch);

    arr_t & H_a = bssn->H_a; // field
    
    const hier::Box& box = bssn->DIFFchi_a_pdata->getGhostBox();
    const hier::Box& inner_box = patch->getBox();

    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];

    const int * inner_lower = &inner_box.lower()[0];
    const int * inner_upper = &inner_box.upper()[0];

    for(int k = lower[2]; k <= upper[2]; k++)
    {
      for(int j = lower[1]; j <= upper[1]; j++)
      {
        for(int i = lower[0]; i <= upper[0]; i++)
        {
          H_a(i, j, k) = sqrt(8.0 * PI * tot_r/3.0);
        }
      }
    }
  }
  return 1;
}

  
bool scalar_ic_set_scalar_gaussian_collapse(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, BSSN * bssn, Scalar * scalar,
  std::shared_ptr<tbox::Database> cosmo_scalar_db)
{
  return 0;
}

  
}
