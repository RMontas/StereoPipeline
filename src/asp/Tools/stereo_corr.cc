// __BEGIN_LICENSE__
//  Copyright (c) 2009-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file stereo_corr.cc
///

#include <vw/InterestPoint.h>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <vw/Stereo/CorrelationView.h>
#include <vw/Stereo/CostFunctions.h>
#include <vw/Stereo/DisparityMap.h>
#include <asp/Tools/stereo.h>
#include <asp/Core/DemDisparity.h>
#include <asp/Core/LocalHomography.h>
#include <asp/Sessions/StereoSession.h>
#include <xercesc/util/PlatformUtils.hpp>

#include <asp/Core/InterestPointMatching.h>
#include <vw/Stereo/StereoModel.h>

#include <asp/Core/AffineEpipolar.cc>
#include <asp/Core/AffineEpipolar.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/InterestPoint/InterestData.h>

using namespace vw;
using namespace vw::stereo;
using namespace asp;
using namespace std;


 /// Homography IP matching // Ricardo Monteiro - return ip matching
  ///
  /// This applies only the homography constraint. Not the best...
  template <class Image1T, class Image2T>
  bool homography_ip_matching1( vw::ImageViewBase<Image1T> const& image1,
			       vw::ImageViewBase<Image2T> const& image2,
			       int ip_per_tile,
			       std::string const& output_name,
			       int inlier_threshold=10,
			       double nodata1 = std::numeric_limits<double>::quiet_NaN(),
			       double nodata2 = std::numeric_limits<double>::quiet_NaN(),
			       std::vector<ip::InterestPoint>& final_ip1 = NULL,
 			       std::vector<ip::InterestPoint>& final_ip2 = NULL );
Vector2i
  affine_epipolar_rectification1( Vector2i const& left_size,
                                 Vector2i const& right_size,
                                 std::vector<ip::InterestPoint> const& ip1,
                                 std::vector<ip::InterestPoint> const& ip2,
                                 Matrix<double>& left_matrix,
                                 Matrix<double>& right_matrix );
Vector2i
  homography_rectification1( bool adjust_left_image_size,
			    Vector2i const& left_size,
			    Vector2i const& right_size,
			    std::vector<ip::InterestPoint> const& left_ip,
			    std::vector<ip::InterestPoint> const& right_ip,
			    vw::Matrix<double>& left_matrix,
			    vw::Matrix<double>& right_matrix,
			    double threshRANSAC, double minAvgDeltaY, BBox2i bbox );
bool check_homography_matrix(Matrix<double>       const& H,
			       std::vector<Vector3> const& left_points,
			       std::vector<Vector3> const& right_points,
			       std::vector<size_t>  const& indices,
			       double minAvgDeltaY, BBox2i bbox);
bool check_homography_matrix(	Matrix<double>       const& left_matrix,
								Matrix<double>       const& right_matrix,
			       				std::vector<Vector3> const& left_points,
			       				std::vector<Vector3> const& right_points,
			       				double minAvgDeltaY, 
								BBox2i bbox);
double calcAverageDeltaY(std::vector<ip::InterestPoint> const& left_points, std::vector<ip::InterestPoint> const& right_points);
double calcAverageDeltaY(std::vector<Vector3> const& left_points, std::vector<Vector3> const& right_points);
BBox2f calcSearchRange(std::vector<ip::InterestPoint> const& left_ip, std::vector<ip::InterestPoint> const& right_ip, Matrix<double> const& left_matrix, Matrix<double> const& right_matrix, double multi);
vw::Matrix<double> piecewiseAlignment(ImageView<float> left_image, 
									ImageView<float> right_image,
									ImageView<float> tile_left_image,
									ImageView<float> tile_right_image,	
									BBox2i bbox);

vw::Matrix<double> piecewiseAlignment_homography(ImageView<float> left_image, 
												ImageView<float> right_image,
												ImageView<float> tile_left_image,
												ImageView<float> tile_right_image,	
												BBox2i bbox);
BBox2f piecewiseAlignment_affineepipolar(	ImageView<float> left_image, 
										ImageView<float> right_image,
										ImageView<float> tile_left_image,
										ImageView<float> tile_right_image,	
										BBox2i bbox,
										Vector2i& left_size,
										Vector2i& right_size,
										vw::Matrix<double>& left_matrix,
										vw::Matrix<double>& right_matrix,
										BBox2f local_search_range);

/// Returns the properly cast cost mode type
stereo::CostFunctionType get_cost_mode_value() {
  switch(stereo_settings().cost_mode) {
    case 0: return stereo::ABSOLUTE_DIFFERENCE;
    case 1: return stereo::SQUARED_DIFFERENCE;
    case 2: return stereo::CROSS_CORRELATION;
    case 3: return stereo::CENSUS_TRANSFORM;
    case 4: return stereo::TERNARY_CENSUS_TRANSFORM;
    default: 
      vw_throw( ArgumentErr() << "Unknown value " << stereo_settings().cost_mode << " for cost-mode.\n" );
  };
}

/// Determine the proper subpixel mode to be used with SGM correlation
SemiGlobalMatcher::SgmSubpixelMode get_sgm_subpixel_mode() {

  switch(stereo_settings().subpixel_mode) {
    case  6: return SemiGlobalMatcher::SUBPIXEL_LINEAR;
    case  7: return SemiGlobalMatcher::SUBPIXEL_POLY4;
    case  8: return SemiGlobalMatcher::SUBPIXEL_COSINE;
    case  9: return SemiGlobalMatcher::SUBPIXEL_PARABOLA;
    case 10: return SemiGlobalMatcher::SUBPIXEL_NONE;
    case 11: return SemiGlobalMatcher::SUBPIXEL_LC_BLEND;
    default: return SemiGlobalMatcher::SUBPIXEL_LC_BLEND;
  };
}


// Read the search range from D_sub, and scale it to the full image
void read_search_range_from_dsub(ASPGlobalOptions & opt){

  // No D_sub is generated or should be used for seed mode 0.
  if (stereo_settings().seed_mode == 0)
    return;

  DiskImageView<vw::uint8> Lmask(opt.out_prefix + "-lMask.tif"),
                           Rmask(opt.out_prefix + "-rMask.tif");

  DiskImageView<PixelGray<float> > left_sub ( opt.out_prefix+"-L_sub.tif" ),
                                   right_sub( opt.out_prefix+"-R_sub.tif" );

  Vector2 downsample_scale( double(left_sub.cols()) / double(Lmask.cols()),
                            double(left_sub.rows()) / double(Lmask.rows()) );

  std::string d_sub_file = opt.out_prefix + "-D_sub.tif";
  if (!fs::exists(d_sub_file))
    return;

  ImageView<PixelMask<Vector2f> > sub_disp;
  read_image(sub_disp, d_sub_file);
  BBox2i search_range = stereo::get_disparity_range( sub_disp );
  search_range.min() = floor(elem_quot(search_range.min(),downsample_scale));
  search_range.max() = ceil (elem_quot(search_range.max(),downsample_scale));
  stereo_settings().search_range = search_range;
  
  vw_out() << "\t--> Read search range from D_sub: " << search_range << "\n";
}



/// Produces the low-resolution disparity file D_sub
void produce_lowres_disparity( ASPGlobalOptions & opt ) {

  // Set up handles to read the input images
  DiskImageView<vw::uint8> Lmask(opt.out_prefix + "-lMask.tif"),
                           Rmask(opt.out_prefix + "-rMask.tif");

  DiskImageView<PixelGray<float> > left_sub ( opt.out_prefix+"-L_sub.tif" ),
                                   right_sub( opt.out_prefix+"-R_sub.tif" );

  DiskImageView<uint8> left_mask_sub ( opt.out_prefix+"-lMask_sub.tif" ),
                       right_mask_sub( opt.out_prefix+"-rMask_sub.tif" );

  Vector2 downsample_scale( double(left_sub.cols()) / double(Lmask.cols()),
                            double(left_sub.rows()) / double(Lmask.rows()) );
  double mean_scale = (downsample_scale[0] + downsample_scale[1]) / 2.0;

  // Compute the initial search range in the subsampled image
  BBox2i search_range( floor(elem_prod(downsample_scale,stereo_settings().search_range.min())),
                       ceil (elem_prod(downsample_scale,stereo_settings().search_range.max())) );

  if ( stereo_settings().seed_mode == 1 ) {

    // Use low-res correlation to get the low-res disparity
    Vector2i expansion( search_range.width(),
                  			search_range.height() );
    expansion *= stereo_settings().seed_percent_pad / 2.0f;
    // Expand by the user selected amount. Default is 25%.
    search_range.min() -= expansion;
    search_range.max() += expansion;
    //VW_OUT(DebugMessage,"asp") << "D_sub search range: " << search_range << " px\n";
    std::cout << "D_sub search range: " << search_range << " px\n";
    stereo::CostFunctionType cost_mode = get_cost_mode_value();
    Vector2i kernel_size  = stereo_settings().corr_kernel;
    int corr_timeout      = 5*stereo_settings().corr_timeout; // 5x, so try hard
    const int rm_half_kernel = 5; // Filter kernel size used by CorrelationView
    double seconds_per_op = 0.0;
    if (corr_timeout > 0)
      seconds_per_op = calc_seconds_per_op(cost_mode, left_sub, right_sub, kernel_size);

    SemiGlobalMatcher::SgmSubpixelMode sgm_subpixel_mode = get_sgm_subpixel_mode();
    Vector2i sgm_search_buffer = stereo_settings().sgm_search_buffer;;

    if (stereo_settings().rm_quantile_multiple <= 0.0)
    {
      // If we can process the entire image in one tile, don't use a collar.
      int collar_size = stereo_settings().sgm_collar_size;
      if ((opt.raster_tile_size[0] > left_sub.cols()) &&
          (opt.raster_tile_size[1] > left_sub.rows())   )
        collar_size = 0;
    
      // Warning: A giant function call approaches!
      // TODO: Why the extra filtering step here? PyramidCorrelationView already performs 1-3 iterations of outlier removal.
      std::string d_sub_file = opt.out_prefix + "-D_sub.tif";
      vw_out() << "Writing: " << d_sub_file << std::endl;
      vw::cartography::block_write_gdal_image( // Write to disk
          d_sub_file,
          rm_outliers_using_thresh( // Throw out individual pixels that are far from any neighbors
              vw::stereo::pyramid_correlate( // Compute image correlation using the PyramidCorrelationView class
                  left_sub, right_sub,
                  left_mask_sub, right_mask_sub,
                  vw::stereo::PREFILTER_LOG, stereo_settings().slogW,
                  search_range, kernel_size, cost_mode,
                  corr_timeout, seconds_per_op,
                  stereo_settings().xcorr_threshold, 
                  stereo_settings().min_xcorr_level,
                  rm_half_kernel,
                  stereo_settings().corr_max_levels,
                  static_cast<vw::stereo::CorrelationAlgorithm>(stereo_settings().stereo_algorithm),
                  collar_size, sgm_subpixel_mode, sgm_search_buffer, stereo_settings().corr_memory_limit_mb,
                  stereo_settings().corr_blob_filter_area*mean_scale,
                  stereo_settings().stereo_debug
              ),
              // To do: all these hard-coded values must be replaced with
              // appropriate params from user's stereo.default, for
              // consistency with how disparity is filtered in stereo_fltr,
              // when invoking disparity_cleanup_using_thresh.
              1, 1, // in stereo.default we have 5 5
              // Changing below the hard-coded value from 2.0 to using a
              // param.  The default value will still be 2.0 but is now
              // modifiable. Need to get rid of the 2.0/3.0 factor and
              // study how it affects the result.
              stereo_settings().rm_threshold*2.0/3.0,
              // Another change of hard-coded value to param. Get rid of 0.5/0.6
              // and study the effect.
              (stereo_settings().rm_min_matches/100.0)*0.5/0.6
          ), // End outlier removal arguments
          opt,
          TerminalProgressCallback("asp", "\t--> Low-resolution disparity:")
      );
      // End of giant function call block
    }
    else { // Use quantile based filtering - This filter needs to be profiled to improve its speed.
    
      // Compute image correlation using the PyramidCorrelationView class
      ImageView< PixelMask<Vector2f> > disp_image = vw::stereo::pyramid_correlate( 
                  left_sub, right_sub,
                  left_mask_sub, right_mask_sub,
                  vw::stereo::PREFILTER_LOG, stereo_settings().slogW,
                  search_range, kernel_size, cost_mode,
                  corr_timeout, seconds_per_op,
                  stereo_settings().xcorr_threshold, 
                  stereo_settings().min_xcorr_level,
                  rm_half_kernel,
                  stereo_settings().corr_max_levels,
                  static_cast<vw::stereo::CorrelationAlgorithm>(stereo_settings().stereo_algorithm), 
                  0, // No collar here, the entire image is written at once.
                  sgm_subpixel_mode, sgm_search_buffer, stereo_settings().corr_memory_limit_mb,
                  0, // Don't combine blob filtering with quantile filtering
                  stereo_settings().stereo_debug
              );

      std::string d_sub_file = opt.out_prefix + "-D_sub.tif";
      vw_out() << "Writing: " << d_sub_file << std::endl;
      vw::cartography::write_gdal_image( // Write to disk while removing outliers
          d_sub_file,
          rm_outliers_using_quantiles( // Throw out individual pixels that are far from any neighbors
              disp_image,
              stereo_settings().rm_quantile_percentile, stereo_settings().rm_quantile_multiple
          ),
          opt,
          TerminalProgressCallback("asp", "\t--> Low-resolution disparity:")
      );
    }

  }else if ( stereo_settings().seed_mode == 2 ) {
    // Use a DEM to get the low-res disparity
    boost::shared_ptr<camera::CameraModel> left_camera_model, right_camera_model;
    opt.session->camera_models(left_camera_model, right_camera_model);
    produce_dem_disparity(opt, left_camera_model, right_camera_model, opt.session->name());
  }else if ( stereo_settings().seed_mode == 3 ) {
    // D_sub is already generated by now by sparse_disp
  }

  read_search_range_from_dsub(opt); // TODO: We already call this when needed!
} // End produce_lowres_disparity


/// Adjust IP lists if alignment matrices are present.
double adjust_ip_for_align_matrix(std::string               const& out_prefix,
                                  vector<ip::InterestPoint>      & ip_left,
                                  vector<ip::InterestPoint>      & ip_right,
                                  double                    const  ip_scale) {

  // Check for alignment files
  bool left_align  = fs::exists(out_prefix+"-align-L.exr");
  bool right_align = fs::exists(out_prefix+"-align-R.exr");
  if (!left_align && !right_align)
    return ip_scale; // No alignment files -> Nothing to do.

  // Load alignment matrices
  Matrix<double> align_left_matrix  = math::identity_matrix<3>();
  Matrix<double> align_right_matrix = math::identity_matrix<3>();
  if (left_align)
    read_matrix(align_left_matrix, out_prefix + "-align-L.exr");
  if (right_align)
    read_matrix(align_right_matrix, out_prefix + "-align-R.exr");

  // Loop through all the IP we found
  for ( size_t i = 0; i < ip_left.size(); i++ ) {
    // Apply the alignment transforms to the recorded IP
    Vector3 l = align_left_matrix  * Vector3(ip_left [i].x, ip_left [i].y, 1);
    Vector3 r = align_right_matrix * Vector3(ip_right[i].x, ip_right[i].y, 1);

    // Normalize the coordinates, but don't divide by 0
    if (l[2] == 0 || r[2] == 0) 
      continue;
    l /= l[2];
    r /= r[2];
    
    ip_left [i].x = l[0];
    ip_left [i].y = l[1];
    ip_right[i].x = r[0];
    ip_right[i].y = r[1];
  }
  return 1.0; // If alignment files are present they take care of the scaling.
	      
} // End adjust_ip_for_align_matrix

  // TODO: Duplicate of hidden function in vw/src/InterestPoint/Matcher.cc!
  std::string strip_path(std::string out_prefix, std::string filename){

    // If filename starts with out_prefix followed by dash, strip both.
    // Also strip filename extension.

    std::string ss = out_prefix + "-";
    size_t found = filename.find(ss);

    if (found != std::string::npos)
      filename.erase(found, ss.length());

    filename = fs::path(filename).stem().string();

    return filename;
  }

/// Detect IP in the _sub images or the original images if they are not too large.
/// - Usually an IP file is written in stereo_pprc, but for some input scenarios
///   this function will need to be used to generate them here.
/// - The input match file path can be changed depending on what exists on disk.
/// - Returns the scale from the image used for IP to the full size image.
/// - The binary interest point file will be written to disk.
double compute_ip(ASPGlobalOptions & opt, std::string & match_filename) {

  vw_out() << "\t    * Loading images for IP detection.\n";

  // Choose whether to use the full or _sub images

  // Use the full image if all dimensions are smaller than this.
  const int SIZE_CUTOFF = 8000;

  const std::string left_image_path_full  = opt.out_prefix+"-L.tif";
  const std::string right_image_path_full = opt.out_prefix+"-R.tif";
  const std::string left_image_path_sub   = opt.out_prefix+"-L_sub.tif";
  const std::string right_image_path_sub  = opt.out_prefix+"-R_sub.tif";

  // TODO: Just call the right function everywhere rather than computing its result by hand.
  const std::string full_match_file       = ip::match_filename(opt.out_prefix, opt.in_file1, opt.in_file2);
  const std::string sub_match_file        = opt.out_prefix + "-L_sub__R_sub.match";
  const std::string aligned_match_file    = opt.out_prefix + "-L__R.match";

  // TODO: The logic below is wrong. Don't read the first match file
  // that happens to exist on disk and hope for the best.  That could
  // be an incorrect one. At this stage we know exactly the files that
  // need processing. Check if the desired file exists, and read that
  // one, or create it if missing.
  
  // Try the full match file first
  if (fs::exists(full_match_file)) {
    vw_out() << "IP file found: " << full_match_file << std::endl;
    match_filename = full_match_file;
    return 1.0;
  }

  // TODO: Unify with function in vw/src/InterestPoint/Matcher.h!
  // filenames longer than this must be chopped, as too long names
  // cause problems later with boost.
  int max_len = 40;
  std::string name1 = strip_path(opt.out_prefix, opt.in_file1).substr(0, max_len);
  std::string name2 = strip_path(opt.out_prefix, opt.in_file2).substr(0, max_len);

  // Next try the cropped match file names which will be at full scale.
  // TODO: This is unnecessary. Just call the right function to find
  // the match file.
  std::vector<std::string> match_names;
  match_names.push_back(opt.out_prefix + "-L-cropped__R-cropped.match");
  match_names.push_back(opt.out_prefix + "-"+name1+"__R-cropped.match");
  match_names.push_back(opt.out_prefix + "-L-cropped__"+name2+".match");
  match_names.push_back(aligned_match_file);
  for (size_t i=0; i<match_names.size(); ++i) {
    if (fs::exists(match_names[i])) {
      vw_out() << "IP file found: " << match_names[i] << std::endl;
      match_filename = match_names[i];
      return 1.0;
    }
  }

  // Now try the sub match file, which requires us to compute the scale.
  std::string left_image_path  = left_image_path_full;
  std::string right_image_path = right_image_path_full;
  Vector2i full_size     = file_image_size(left_image_path_full);
  bool     use_full_size = (((full_size[0] < SIZE_CUTOFF) && (full_size[1] < SIZE_CUTOFF))
                            || ((stereo_settings().alignment_method != "epipolar") &&
                                (stereo_settings().alignment_method != "none"    )   ));
  // Other alignment methods find IP in the stereo_pprc phase using the full size.

  // Compute the scale.
  double ip_scale = 1.0;
  if (!use_full_size) {
    left_image_path  = left_image_path_sub;
    right_image_path = right_image_path_sub;
    
    ip_scale = sum(elem_quot( Vector2(file_image_size( opt.out_prefix+"-L_sub.tif" )),
                              Vector2(file_image_size( opt.out_prefix+"-L.tif" ) ) )) +
               sum(elem_quot( Vector2(file_image_size( opt.out_prefix+"-R_sub.tif" )),
                              Vector2(file_image_size( opt.out_prefix+"-R.tif" ) ) ));
    ip_scale /= 4.0f;
    match_filename = sub_match_file; // If not using full size we should expect this file
    
    // Check for the file.
    if (fs::exists(sub_match_file)) {
      vw_out() << "IP file found: " << sub_match_file << std::endl;
      return ip_scale;
    }
  }
  else
    match_filename = aligned_match_file;

  vw_out() << "No IP file found, computing IP now.\n";
  
  // Load the images
  boost::shared_ptr<DiskImageResource> left_rsrc (DiskImageResourcePtr(left_image_path )),
                                       right_rsrc(DiskImageResourcePtr(right_image_path));

  // Read the no-data values written to disk previously when
  // the normalized left and right sub-images were created.
  float left_nodata_value  = numeric_limits<float>::quiet_NaN();
  float right_nodata_value = numeric_limits<float>::quiet_NaN();
  if (left_rsrc->has_nodata_read ()) left_nodata_value  = left_rsrc->nodata_read();
  if (right_rsrc->has_nodata_read()) right_nodata_value = right_rsrc->nodata_read();
  
  // These images should be small enough to fit in memory
  ImageView<float> left_image  = DiskImageView<float>(left_rsrc);
  ImageView<float> right_image = DiskImageView<float>(right_rsrc);

  // No interest point operations have been performed before
  vw_out() << "\t    * Locating Interest Points\n";

  // Use this code in a relatively specific case
  // - Only tested with IceBridge data so far!
  // - Some changes will be required for this to work in more general cases.
  bool success = false;
  if (use_full_size && opt.session->is_nadir_facing() && 
      (stereo_settings().alignment_method == "epipolar") ) {

    // Load camera models
    boost::shared_ptr<camera::CameraModel> left_camera_model, right_camera_model;
    opt.session->camera_models(left_camera_model, right_camera_model);
    
    // Obtain the datum
    const bool use_sphere_for_isis = false;
    cartography::Datum datum = opt.session->get_datum(left_camera_model.get(), use_sphere_for_isis);

    // Since these are epipolar aligned images it should be small
    double epipolar_threshold = 5;
    if (stereo_settings().epipolar_threshold > 0)
      epipolar_threshold = stereo_settings().epipolar_threshold;

    const bool single_threaded_camera = false;
    success = ip_matching(single_threaded_camera,
                          left_camera_model.get(), right_camera_model.get(),
                          left_image, right_image,
                          stereo_settings().ip_per_tile,
                          datum, match_filename, epipolar_threshold,
                          stereo_settings().ip_uniqueness_thresh,
                          left_nodata_value, right_nodata_value);
  } // End nadir epipolar full image case
  else {
    // In all other cases, run a more general IP matcher.
    
    // TODO: Depending on alignment method, we can tailor the IP filtering strategy.
    double thresh_factor = stereo_settings().ip_inlier_factor; // 1/15 by default

    // This range is extra large to handle elevation differences.
    const int inlier_threshold = 200*(15.0*thresh_factor);  // 200 by default
    
    success = asp::homography_ip_matching(left_image, right_image,
                                          stereo_settings().ip_per_tile,
                                          match_filename, inlier_threshold,
                                          left_nodata_value, right_nodata_value);
  }

  if (!success)
    vw_throw(ArgumentErr() << "Could not find interest points.\n");

  return ip_scale;
}




// TODO: Move this histogram code!  Merge with image histogram code!


  /// Compute simple histogram from a vector of data
  void histogram( std::vector<double> const& values, int num_bins, double min_val, double max_val,
                  std::vector<double> &hist, std::vector<double> &bin_centers){
    
    VW_ASSERT(num_bins > 0, ArgumentErr() << "histogram: number of input bins must be positive");
    
    // TODO: Verify max/min values!
    
    // Populate the list of bin centers
    // The min and max vals represent the outer limits of the available bins.
    const double range     = max_val - min_val;
    const double bin_width = range / num_bins;
    bin_centers.resize(num_bins);
    for (int i=0; i<num_bins; ++i)
      bin_centers[i] = min_val + i*bin_width + bin_width/2.0;
    
    hist.assign(num_bins, 0.0);
    for (size_t i = 0; i < values.size(); ++i){
      double val = values[i];
      int bin = (int)round( (num_bins - 1) * ((val - min_val)/range) );
      
      // Saturate the bin assignment to prevent a memory exception.
      if (bin < 0)
        bin = 0;
      if (bin > num_bins-1)
        bin = num_bins-1;
        
      ++(hist[bin]);
    }
    return;
  }

/// Use existing interest points to compute a search range
/// - This function could use improvement!
/// - Should it be used in all cases?
BBox2i approximate_search_range(ASPGlobalOptions & opt, 
                                double ip_scale, std::string const& match_filename) {

  vw_out() << "\t--> Using interest points to determine search window.\n";
  vector<ip::InterestPoint> in_ip1, in_ip2, matched_ip1, matched_ip2;

  // The interest points must have been created outside this function
  if (!fs::exists(match_filename))
    vw_throw( ArgumentErr() << "Missing IP file: " << match_filename);

  vw_out() << "\t    * Loading match file: " << match_filename << "\n";
  ip::read_binary_match_file(match_filename, in_ip1, in_ip2);

  // Handle alignment matrices if they are present
  // - Scale is reset to 1.0 if alignment matrices are present.
  ip_scale = adjust_ip_for_align_matrix(opt.out_prefix, in_ip1, in_ip2, ip_scale);
  vw_out() << "\t    * IP computed at scale: " << ip_scale << "\n";
  float i_scale = 1.0/ip_scale;

  // Filter out IPs which fall outside the specified elevation range
  boost::shared_ptr<camera::CameraModel> left_camera_model, right_camera_model;
  opt.session->camera_models(left_camera_model, right_camera_model);
  cartography::Datum datum = opt.session->get_datum(left_camera_model.get(), false);

  // We already corrected for align matrix, so transforms should be identity here.
  vw::TransformRef left_tx  = vw::TransformRef(vw::TranslateTransform(0,0));
  vw::TransformRef right_tx = vw::TransformRef(vw::TranslateTransform(0,0));

  // Filter out IPs which fall outside the specified elevation and lonlat range
  // TODO: Don't do this with cropped input images!!!!!
  size_t num_left = asp::filter_ip_by_lonlat_and_elevation(left_camera_model.get(),
							   right_camera_model.get(),
							   datum, in_ip1, in_ip2,
							   left_tx, right_tx, ip_scale,
							   stereo_settings().elevation_limit,
							   stereo_settings().lon_lat_limit,
							   matched_ip1, matched_ip2);


  if (num_left == 0)
    vw_throw(ArgumentErr() << "No IPs left after elevation filtering!");

  // Find search window based on interest point matches

  // Record the disparities for each point pair
  size_t num_ip = matched_ip1.size();
  std::vector<double> dx, dy;
  double min_dx, min_dy, max_dx, max_dy;
  min_dx = min_dy = std::numeric_limits<double>::max();
  max_dx = max_dy = std::numeric_limits<double>::min();
  for (size_t i = 0; i < num_ip; i++) {
    double diffX = i_scale * (matched_ip2[i].x - matched_ip1[i].x);
    double diffY = i_scale * (matched_ip2[i].y - matched_ip1[i].y);      
    dx.push_back(diffX);
    dy.push_back(diffY);
    
    if (diffX < min_dx) min_dx = diffX;  // Could be smarter about how this is done.
    if (diffY < min_dy) min_dy = diffY;
    if (diffX > max_dx) max_dx = diffX;
    if (diffY > max_dy) max_dy = diffY;
  }
  num_ip = dx.size();
  
  // Compute histograms
  const int NUM_BINS = 2000; // Accuracy is important with scaled pixels
  std::vector<double> hist_x, centers_x, hist_y, centers_y;
  histogram(dx, NUM_BINS, min_dx, max_dx, hist_x, centers_x);
  histogram(dy, NUM_BINS, min_dy, max_dy, hist_y, centers_y);
  
  //printf("min x,y = %lf, %lf, max x,y = %lf, %lf\n", min_dx, min_dy, max_dx, max_dy);
  
  // Compute search ranges
  const double MAX_PERCENTILE = 0.95;
  const double MIN_PERCENTILE = 0.05;
  double search_scale = 2.0;
  size_t min_bin_x = get_histogram_percentile(hist_x, MIN_PERCENTILE);
  size_t min_bin_y = get_histogram_percentile(hist_y, MIN_PERCENTILE);
  size_t max_bin_x = get_histogram_percentile(hist_x, MAX_PERCENTILE);
  size_t max_bin_y = get_histogram_percentile(hist_y, MAX_PERCENTILE);
  Vector2 search_min(centers_x[min_bin_x],
                     centers_y[min_bin_y]);
  Vector2 search_max(centers_x[max_bin_x],
                     centers_y[max_bin_y]);
  //std::cout << "Unscaled search range = " << BBox2i(search_min, search_max) << std::endl;
  Vector2 search_center = (search_max + search_min) / 2.0;
  //std::cout << "search_center = " << search_center << std::endl;
  Vector2 d_min = search_min - search_center; // TODO: Make into a bbox function!
  Vector2 d_max = search_max - search_center;
  //std::cout << "d_min = " << d_min << std::endl;
  //std::cout << "d_max = " << d_max << std::endl;
  search_min = d_min*search_scale + search_center;
  search_max = d_max*search_scale + search_center;
  //std::cout << "Scaled search range = " << BBox2i(search_min, search_max) << std::endl;
  
/*
   // Debug code to print all the points
  for (size_t i = 0; i < matched_ip1.size(); i++) {
    Vector2f diff(i_scale * (matched_ip2[i].x - matched_ip1[i].x), 
                  i_scale * (matched_ip2[i].y - matched_ip1[i].y));
    //Vector2f diff(matched_ip2[i].x - matched_ip1[i].x, 
    //              matched_ip2[i].y - matched_ip1[i].y);
    vw_out(InfoMessage,"asp") << matched_ip1[i].x <<", "<<matched_ip1[i].y 
              << " <> " 
              << matched_ip2[i].x <<", "<<matched_ip2[i].y 
               << " DIFF " << diff << " DIFF-M " << diff-mean << endl;
  }
*/
  
  //vw_out(InfoMessage,"asp") << "i_scale is : "       << i_scale << endl;
  BBox2i search_range(search_min, search_max);
  return search_range;
} // End function approximate_search_range


/// The first step of correlation computation.
void lowres_correlation( ASPGlobalOptions & opt ) {

  vw_out() << "\n[ " << current_posix_time_string() << " ] : Stage 1 --> LOW-RESOLUTION CORRELATION \n";

  // Working out search range if need be
  if (stereo_settings().is_search_defined()) {
    vw_out() << "\t--> Using user-defined search range.\n";

    // Update user provided search range based on input crops
    bool crop_left  = (stereo_settings().left_image_crop_win  != BBox2i(0, 0, 0, 0));
    bool crop_right = (stereo_settings().right_image_crop_win != BBox2i(0, 0, 0, 0));
    if (crop_left && !crop_right)
      stereo_settings().search_range += stereo_settings().left_image_crop_win.min();
    if (!crop_left && crop_right)
      stereo_settings().search_range -= stereo_settings().right_image_crop_win.min();

  }else if (stereo_settings().seed_mode == 2){
    // Do nothing as we will compute the search range based on D_sub
  }else if (stereo_settings().seed_mode == 3){
    // Do nothing as low-res disparity (D_sub) is already provided by sparse_disp
  } else { // Regular seed mode

    // If there is no match file for the input images, gather some IP from the
    // low resolution images. This routine should only run for:
    //   Pinhole + Epipolar
    //   Alignment method none
    //   Cases where either input image is cropped, in which case the IP name is different.
    // Everything else should gather IP's all the time during stereo_pprc.
    // - TODO: When inputs are cropped, use the cropped IP!
    

    // Compute new IP and write them to disk.
    // - If IP are already on disk this function will load them instead.
    // - This function will choose an appropriate IP computation based on the input images.
    string match_filename;
    double ip_scale;
    ip_scale = compute_ip(opt, match_filename);

    // This function applies filtering to find good points
    stereo_settings().search_range = approximate_search_range(opt, ip_scale, match_filename);
  
 
    vw_out() << "\t--> Detected search range: " << stereo_settings().search_range << "\n";
  } // End of case where we had to calculate the search range

  // If the user specified a search range limit, apply it here.
  if ((stereo_settings().search_range_limit.min() != Vector2i()) || 
      (stereo_settings().search_range_limit.max() != Vector2i())   ) {     
    stereo_settings().search_range.crop(stereo_settings().search_range_limit);
    vw_out() << "\t--> Detected search range constrained to: " << stereo_settings().search_range << "\n";
  }

  // At this point stereo_settings().search_range is populated

  DiskImageView<vw::uint8> Lmask(opt.out_prefix + "-lMask.tif"),
                           Rmask(opt.out_prefix + "-rMask.tif");

  // Performing disparity on sub images
  if ( stereo_settings().seed_mode > 0 ) {

    // Reuse prior existing D_sub if it exists, unless we
    // are cropping the images each time, when D_sub must
    // be computed anew each time.
    bool crop_left  = (stereo_settings().left_image_crop_win  != BBox2i(0, 0, 0, 0));
    bool crop_right = (stereo_settings().right_image_crop_win != BBox2i(0, 0, 0, 0));
    bool rebuild    = crop_left || crop_right;

    string sub_disp_file = opt.out_prefix+"-D_sub.tif";
    try {
      vw_log().console_log().rule_set().add_rule(-1,"fileio");
      DiskImageView<PixelMask<Vector2f> > test(sub_disp_file);
      vw_settings().reload_config();
    } catch (vw::IOErr const& e) {
      vw_settings().reload_config();
      rebuild = true;
    } catch (vw::ArgumentErr const& e ) {
      // Throws on a corrupted file.
      vw_settings().reload_config();
      rebuild = true;
    }

    if ( rebuild )
      produce_lowres_disparity(opt); // Note: This does not always remake D_sub!
    else
      vw_out() << "\t--> Using cached low-resolution disparity: " << sub_disp_file << "\n";
  }

  // Create the local homographies based on D_sub
  if (stereo_settings().seed_mode > 0 && stereo_settings().use_local_homography){
    string local_hom_file = opt.out_prefix + "-local_hom.txt";
    try {
      ImageView<Matrix3x3> local_hom;
      read_local_homographies(local_hom_file, local_hom);
    } catch (vw::IOErr const& e) {
      create_local_homographies(opt);
    }
  }

  vw_out() << "\n[ " << current_posix_time_string() << " ] : LOW-RESOLUTION CORRELATION FINISHED \n";
} // End lowres_correlation



/// This correlator takes a low resolution disparity image as an input
/// so that it may narrow its search range for each tile that is processed.
class SeededCorrelatorView : public ImageViewBase<SeededCorrelatorView> {
  DiskImageView<PixelGray<float> >   m_left_image;
  DiskImageView<PixelGray<float> >   m_right_image;
  DiskImageView<vw::uint8> m_left_mask;
  DiskImageView<vw::uint8> m_right_mask;
  ImageViewRef<PixelMask<Vector2f> > m_sub_disp;
  ImageViewRef<PixelMask<Vector2i> > m_sub_disp_spread;
  ImageView<Matrix3x3> & m_local_hom;

  // Settings
  Vector2  m_upscale_factor;
  BBox2i   m_seed_bbox;
  Vector2i m_kernel_size;
  stereo::CostFunctionType m_cost_mode;
  int      m_corr_timeout;
  double   m_seconds_per_op;

public:

  // Set these input types here instead of making them template arguments
  typedef DiskImageView<PixelGray<float> >   ImageType;
  typedef DiskImageView<vw::uint8>           MaskType;
  typedef ImageViewRef<PixelMask<Vector2f> > DispSeedImageType;
  typedef ImageViewRef<PixelMask<Vector2i> > SpreadImageType;
  typedef ImageType::pixel_type InputPixelType;

  SeededCorrelatorView( ImageType             const& left_image,
                        ImageType             const& right_image,
                        MaskType              const& left_mask,
                        MaskType              const& right_mask,
                        DispSeedImageType     const& sub_disp,
                        SpreadImageType       const& sub_disp_spread,
                        ImageView<Matrix3x3>  & local_hom,
                        Vector2i const& kernel_size,
                        stereo::CostFunctionType cost_mode,
                        int corr_timeout, double seconds_per_op):
    m_left_image(left_image.impl()), m_right_image(right_image.impl()),
    m_left_mask (left_mask.impl ()), m_right_mask (right_mask.impl ()),
    m_sub_disp(sub_disp.impl()), m_sub_disp_spread(sub_disp_spread.impl()),
    m_local_hom(local_hom),
    m_kernel_size(kernel_size),  m_cost_mode(cost_mode),
    m_corr_timeout(corr_timeout), m_seconds_per_op(seconds_per_op){ 
    m_upscale_factor[0] = double(m_left_image.cols()) / m_sub_disp.cols();
    m_upscale_factor[1] = double(m_left_image.rows()) / m_sub_disp.rows();
    m_seed_bbox = bounding_box( m_sub_disp );
  }

  // Image View interface
  typedef PixelMask<Vector2f> pixel_type;
  typedef pixel_type          result_type;
  typedef ProceduralPixelAccessor<SeededCorrelatorView> pixel_accessor;

  inline int32 cols  () const { return m_left_image.cols(); }
  inline int32 rows  () const { return m_left_image.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  inline pixel_type operator()(double /*i*/, double /*j*/, int32 /*p*/ = 0) const {
    vw_throw(NoImplErr() << "SeededCorrelatorView::operator()(...) is not implemented");
    return pixel_type();
  }

  /// Does the work
  typedef CropView<ImageView<pixel_type> > prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {
cout << "start of tile " << bbox << endl;
    bool use_local_homography = stereo_settings().use_local_homography;
    Matrix<double> lowres_hom  = math::identity_matrix<3>();
    Matrix<double> fullres_hom = math::identity_matrix<3>();
    ImageViewRef<InputPixelType> right_trans_img;
    ImageViewRef<vw::uint8     > right_trans_mask;
// piecewise alignment - Ricardo Monteiro
	char outputName[30];
	int ts = ASPGlobalOptions::corr_tile_size();
	int W = bbox.min().x()/ts;
	int H = bbox.min().y()/ts;
	cartography::GdalWriteOptions geo_opt;
    ImageViewRef<InputPixelType> left_trans_img;
    ImageViewRef<vw::uint8     > left_trans_mask;
	int margin = 50;
	BBox2i newBBox = BBox2i(bbox.min().x(), bbox.min().y(), bbox.max().x(), bbox.max().y());
	newBBox.expand(margin);
	newBBox.crop(bounding_box(m_left_image));
	ImageView<PixelGray<float> > tile_right_image = crop(m_right_image.impl(), newBBox);
	ImageView<PixelGray<float> > tile_left_image = crop(m_left_image.impl(), newBBox);
	ImageView<vw::uint8> tile_right_image_mask = crop(m_right_mask.impl(), newBBox);
	ImageView<vw::uint8> tile_left_image_mask = crop(m_left_mask.impl(), newBBox);
	Matrix<double>  align_left_matrix  = math::identity_matrix<3>(),
                   	align_right_matrix = math::identity_matrix<3>();

    bool do_round = true; // round integer disparities after transform

    // User strategies
    BBox2f local_search_range;
    if ( stereo_settings().seed_mode > 0 ) {

      // The low-res version of bbox
      BBox2i seed_bbox( elem_quot(bbox.min(), m_upscale_factor),
			                  elem_quot(bbox.max(), m_upscale_factor) );
      seed_bbox.expand(1);
      seed_bbox.crop( m_seed_bbox );
      // Get the disparity range in d_sub corresponding to this tile.
      VW_OUT(DebugMessage, "stereo") << "Getting disparity range for : " << seed_bbox << "\n";
      DispSeedImageType disparity_in_box = crop( m_sub_disp, seed_bbox );

      if (!use_local_homography){
        local_search_range = stereo::get_disparity_range( disparity_in_box );
      }else{ // use local homography
        int ts = ASPGlobalOptions::corr_tile_size();
        lowres_hom = m_local_hom(bbox.min().x()/ts, bbox.min().y()/ts);
        local_search_range = stereo::get_disparity_range
          (transform_disparities(do_round, seed_bbox,
			     lowres_hom, disparity_in_box)); 
      }

      bool has_sub_disp_spread = ( m_sub_disp_spread.cols() != 0 &&
			                             m_sub_disp_spread.rows() != 0 );
      // Sanity check: If m_sub_disp_spread was provided, it better have the same size as sub_disp.
      if ( has_sub_disp_spread &&
           m_sub_disp_spread.cols() != m_sub_disp.cols() &&
           m_sub_disp_spread.rows() != m_sub_disp.rows() ){
        vw_throw( ArgumentErr() << "stereo_corr: D_sub and D_sub_spread must have equal sizes.\n");
      }

      if (has_sub_disp_spread){
        // Expand the disparity range by m_sub_disp_spread.
        SpreadImageType spread_in_box = crop( m_sub_disp_spread, seed_bbox );

        if (!use_local_homography){
          BBox2f spread = stereo::get_disparity_range( spread_in_box );
          local_search_range.min() -= spread.max();
          local_search_range.max() += spread.max();
        }else{
          DispSeedImageType upper_disp = transform_disparities(do_round, seed_bbox, lowres_hom,
                                                               disparity_in_box + spread_in_box);
          DispSeedImageType lower_disp = transform_disparities(do_round, seed_bbox, lowres_hom,
                                                               disparity_in_box - spread_in_box);
          BBox2f upper_range = stereo::get_disparity_range(upper_disp);
          BBox2f lower_range = stereo::get_disparity_range(lower_disp);

          local_search_range = upper_range;
          local_search_range.grow(lower_range);
        } //endif use_local_homography
      } //endif has_sub_disp_spread

      if (use_local_homography){
        Vector3 upscale(     m_upscale_factor[0],     m_upscale_factor[1], 1 );
        Vector3 dnscale( 1.0/m_upscale_factor[0], 1.0/m_upscale_factor[1], 1 );
        fullres_hom = diagonal_matrix(upscale)*lowres_hom*diagonal_matrix(dnscale);
////////// Ricardo Monteiro - code to overwrite fullres_hom
		local_search_range = stereo::get_disparity_range( disparity_in_box );		
		sprintf(outputName, "tile_R_%d_%d.tif", H, W);
		block_write_gdal_image(outputName, tile_right_image, geo_opt);
		sprintf(outputName, "tile_L_%d_%d.tif", H, W);
		block_write_gdal_image(outputName, tile_left_image, geo_opt);
		Vector2i left_size = newBBox.size();
 		Vector2i right_size = newBBox.size();
		cout << "[tile(" << H << "," << W << " left_size = " << left_size << endl;
		cout << "[tile(" << H << "," << W << " right_size = " << right_size << endl;
	    local_search_range = piecewiseAlignment_affineepipolar(m_left_image.impl(), m_right_image.impl(), tile_left_image.impl(), tile_right_image.impl(), newBBox, left_size, right_size, align_left_matrix, align_right_matrix, local_search_range);
		cout << "[tile(" << H << "," << W << " local_search_range after piecewise alignment = " << local_search_range << endl;
		right_size = left_size;
		cout << "[tile(" << H << "," << W << " left_size after piecewise alignment = " << left_size << endl;
		cout << "[tile(" << H << "," << W << " right_size after piecewise alignment = " << right_size << endl;
		fullres_hom = align_right_matrix;
		m_local_hom(bbox.min().x()/ts, bbox.min().y()/ts) = fullres_hom;
		cout << "[tile(" << H << "," << W << " local_search_range = " << local_search_range << endl;
		
		cout << "[tile(" << H << "," << W << ") " << fullres_hom << "]" << endl;
		cout << "[tile(" << H << "," << W << ") " << align_left_matrix << "]" << endl;
		ImageViewRef< PixelMask<InputPixelType> >
          left_trans_masked_img
          //= transform (copy_mask( m_left_image.impl(),
			//          create_mask(m_left_mask.impl()) ),
			= transform (copy_mask( tile_left_image.impl(),
						create_mask(tile_left_image_mask.impl()) ),
	               HomographyTransform(align_left_matrix),
	              // m_left_image.impl().cols(), m_left_image.impl().rows());
					left_size.x(), left_size.y()); 
        left_trans_img  = apply_mask(left_trans_masked_img);
        left_trans_mask = channel_cast_rescale<uint8>(select_channel(left_trans_masked_img, 1));
/////
        ImageViewRef< PixelMask<InputPixelType> >
          right_trans_masked_img
          //= transform (copy_mask( m_right_image.impl(),
		//	          create_mask(m_right_mask.impl()) ),
		  = transform (copy_mask(tile_right_image.impl(),
					 create_mask(tile_right_image_mask.impl()) ),
	               HomographyTransform(fullres_hom),
	               //m_left_image.impl().cols(), m_left_image.impl().rows());
					right_size.x(), right_size.y()); 
        right_trans_img  = apply_mask(right_trans_masked_img);
        right_trans_mask = channel_cast_rescale<uint8>(select_channel(right_trans_masked_img, 1));

//// write ind tiles
	sprintf(outputName, "piecewiseHomography_R_%d_%d.tif", H, W);
	block_write_gdal_image(outputName, right_trans_img, geo_opt);
	sprintf(outputName, "piecewiseHomography_L_%d_%d.tif", H, W);
	block_write_gdal_image(outputName, left_trans_img, geo_opt);
/////

      } //endif use_local_homography

      local_search_range = grow_bbox_to_int(local_search_range);
      // Expand local_search_range by 1. This is necessary since
      // m_sub_disp is integer-valued, and perhaps the search
      // range was supposed to be a fraction of integer bigger.
      local_search_range.expand(1);
      
      // Scale the search range to full-resolution
      local_search_range.min() = floor(elem_prod(local_search_range.min(),m_upscale_factor));
      local_search_range.max() = ceil (elem_prod(local_search_range.max(),m_upscale_factor));

      // If the user specified a search range limit, apply it here.
      if ((stereo_settings().search_range_limit.min() != Vector2i()) || 
          (stereo_settings().search_range_limit.max() != Vector2i())   ) {     
        local_search_range.crop(stereo_settings().search_range_limit);
        vw_out() << "\t--> Local search range constrained to: " << local_search_range << "\n";
      }

      VW_OUT(DebugMessage, "stereo") << "SeededCorrelatorView("
				     << bbox << ") local search range "
				     << local_search_range << " vs "
				     << stereo_settings().search_range << "\n";

    } else{ // seed mode == 0
      local_search_range = stereo_settings().search_range;
      VW_OUT(DebugMessage,"stereo") << "Searching with " << stereo_settings().search_range << "\n";
    }

    SemiGlobalMatcher::SgmSubpixelMode sgm_subpixel_mode = get_sgm_subpixel_mode();
    Vector2i sgm_search_buffer = stereo_settings().sgm_search_buffer;

    // Now we are ready to actually perform correlation
    const int rm_half_kernel = 5; // Filter kernel size used by CorrelationView
    if (use_local_homography){
      //typedef vw::stereo::PyramidCorrelationView<ImageType, ImageViewRef<InputPixelType>, 
      //                                           MaskType,  ImageViewRef<vw::uint8     > > CorrView;
    //  CorrView corr_view( m_left_image,   right_trans_img,
    //                      m_left_mask,    right_trans_mask,
		typedef vw::stereo::PyramidCorrelationView<ImageViewRef<InputPixelType>, ImageViewRef<InputPixelType>, 
                                                   ImageViewRef<vw::uint8     >, ImageViewRef<vw::uint8     > > CorrView;
	  CorrView corr_view( left_trans_img,   right_trans_img,
                          left_trans_mask,    right_trans_mask,
                          static_cast<vw::stereo::PrefilterModeType>(stereo_settings().pre_filter_mode),
                          stereo_settings().slogW,
                          local_search_range,
                          m_kernel_size,  m_cost_mode,
                          m_corr_timeout, m_seconds_per_op,
                          stereo_settings().xcorr_threshold,
                          stereo_settings().min_xcorr_level,
                          rm_half_kernel,
                          stereo_settings().corr_max_levels,
                          static_cast<vw::stereo::CorrelationAlgorithm>(stereo_settings().stereo_algorithm), 
                          stereo_settings().sgm_collar_size,
                          sgm_subpixel_mode, sgm_search_buffer, stereo_settings().corr_memory_limit_mb,
                          stereo_settings().corr_blob_filter_area,
                          stereo_settings().stereo_debug );
      cout << "end of tile " << newBBox << endl;
      //return corr_view.prerasterize(bbox);
      ImageView<pixel_type> stereo_result = corr_view.prerasterize(bounding_box(left_trans_img));
      ImageView<pixel_type> stereo_result_inv;
      ImageView<vw::uint8     > stereo_result_mask_inv;  
	  ImageView<vw::uint8     > stereo_result_mask = left_trans_mask;
	  // write stereo result
	  sprintf(outputName, "stereo_%d_%d.tif", H, W);
	  block_write_gdal_image(outputName, stereo_result, geo_opt);
	  
		ImageView< PixelMask<pixel_type> >
          stereo_result_masked_img_inv
			= transform (copy_mask(stereo_result.impl(), stereo_result_mask.impl()),
	               HomographyTransform(inverse(align_left_matrix)),
					//tile_left_image.cols(), tile_left_image.rows()); 
					newBBox.width(), newBBox.height());
        stereo_result_inv  = apply_mask(stereo_result_masked_img_inv);
        stereo_result_mask_inv = channel_cast_rescale<uint8>(select_channel(stereo_result_masked_img_inv, 2));
	sprintf(outputName, "stereoINV_%d_%d.tif", H, W);
	//ImageView<pixel_type> image(newBBox.width(), newBBox.height());	
	block_write_gdal_image(outputName, stereo_result_inv, geo_opt);  

	ImageView<pixel_type> stereo_result_corrected(bbox.width(), bbox.height());
	for(int j=0; j<bbox.height(); j++ ){
		for(int i=0; i<bbox.width(); i++ ){
			Vector2 pixel_L_prime = HomographyTransform(align_left_matrix).forward(Vector2(i,j));
			//cout << "pixel_L_prime = " << pixel_L_prime << endl;
			float dx = stereo_result_inv(i+margin,j+margin)[0];
			float dy = stereo_result_inv(i+margin,j+margin)[1];
			Vector2 pixel_R_prime = pixel_L_prime + Vector2(dx,dy);
			//cout << "pixel_R_prime = " << pixel_R_prime << endl;
			Vector2 new_disp = pixel_R_prime - Vector2(i,j);
			//cout << "new_disp = " << new_disp << endl;
			//cout << i-marginMinX << " " << j-marginMinY << " ..  ";
			stereo_result_corrected(i,j)[0] = new_disp.x();
			stereo_result_corrected(i,j)[1] = new_disp.y();
			if(stereo_result_mask_inv(i+margin,j+margin))
				validate(stereo_result_corrected(i,j));
				//stereo_result_corrected(i,j)[2] = 2147483647;
		}
	}

	  //return prerasterize_type(image,-bbox.min().x(),-bbox.min().y(),cols(),rows() );
	return prerasterize_type(stereo_result_corrected,-bbox.min().x(),-bbox.min().y(),cols(),rows() );
    }else{
      typedef vw::stereo::PyramidCorrelationView<ImageType, ImageType, MaskType, MaskType > CorrView;
      CorrView corr_view( m_left_image,   m_right_image,
                          m_left_mask,    m_right_mask,
                          static_cast<vw::stereo::PrefilterModeType>(stereo_settings().pre_filter_mode),
                          stereo_settings().slogW,
                          local_search_range,
                          m_kernel_size,  m_cost_mode,
                          m_corr_timeout, m_seconds_per_op,
                          stereo_settings().xcorr_threshold,
                          stereo_settings().min_xcorr_level,
                          rm_half_kernel,
                          stereo_settings().corr_max_levels,
                          static_cast<vw::stereo::CorrelationAlgorithm>(stereo_settings().stereo_algorithm), 
                          stereo_settings().sgm_collar_size,
                          sgm_subpixel_mode, sgm_search_buffer, stereo_settings().corr_memory_limit_mb,
                          stereo_settings().corr_blob_filter_area,
                          stereo_settings().stereo_debug );
      return corr_view.prerasterize(bbox);
    }
    
  } // End function prerasterize_helper

  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
}; // End class SeededCorrelatorView


/// Main stereo correlation function, called after parsing input arguments.
void stereo_correlation( ASPGlobalOptions& opt ) {

  // The first thing we will do is compute the low-resolution correlation.

  // Note that even when we are told to skip low-resolution correlation,
  // we must still go through the motions when seed_mode is 0, to be
  // able to get a search range, even though we don't write D_sub then.
  if (!stereo_settings().skip_low_res_disparity_comp || stereo_settings().seed_mode == 0)
    lowres_correlation(opt);

  if (stereo_settings().compute_low_res_disparity_only) 
    return; // Just computed the low-res disparity, so quit.

  vw_out() << "\n[ " << current_posix_time_string() << " ] : Stage 1 --> CORRELATION \n";

  read_search_range_from_dsub(opt);

  // If the user specified a search range limit, apply it here.
  if ((stereo_settings().search_range_limit.min() != Vector2i()) || 
      (stereo_settings().search_range_limit.max() != Vector2i())   ) {     
    stereo_settings().search_range.crop(stereo_settings().search_range_limit);
    vw_out() << "\t--> Detected search range constrained to: " << stereo_settings().search_range << "\n";
  }


  // Provide the user with some feedback of what we are actually going to use.
  vw_out()   << "\t--------------------------------------------------\n";
  vw_out()   << "\t   Kernel Size:    " << stereo_settings().corr_kernel << endl;
  if ( stereo_settings().seed_mode > 0 )
    vw_out() << "\t   Refined Search: " << stereo_settings().search_range << endl;
  else
    vw_out() << "\t   Search Range:   " << stereo_settings().search_range << endl;
  vw_out()   << "\t   Cost Mode:      " << stereo_settings().cost_mode << endl;
  vw_out(DebugMessage) << "\t   XCorr Threshold: " << stereo_settings().xcorr_threshold << endl;
  vw_out(DebugMessage) << "\t   Prefilter:       " << stereo_settings().pre_filter_mode << endl;
  vw_out(DebugMessage) << "\t   Prefilter Size:  " << stereo_settings().slogW << endl;
  vw_out() << "\t--------------------------------------------------\n";

  // Load up for the actual native resolution processing
  DiskImageView<PixelGray<float> > left_disk_image (opt.out_prefix+"-L.tif"),
                                   right_disk_image(opt.out_prefix+"-R.tif");
  DiskImageView<vw::uint8> Lmask(opt.out_prefix + "-lMask.tif"),
                           Rmask(opt.out_prefix + "-rMask.tif");
  ImageViewRef<PixelMask<Vector2f> > sub_disp;
  std::string dsub_file   = opt.out_prefix+"-D_sub.tif";
  std::string spread_file = opt.out_prefix+"-D_sub_spread.tif";
  
  if ( stereo_settings().seed_mode > 0 )
    sub_disp = DiskImageView<PixelMask<Vector2f> >(dsub_file);
  ImageViewRef<PixelMask<Vector2i> > sub_disp_spread;
  if ( stereo_settings().seed_mode == 2 ||  stereo_settings().seed_mode == 3 ){
    // D_sub_spread is mandatory for seed_mode 2 and 3.
    sub_disp_spread = DiskImageView<PixelMask<Vector2i> >(spread_file);
  }else if ( stereo_settings().seed_mode == 1 ){
    // D_sub_spread is optional for seed_mode 1, we use it only if it is provided.
    if (fs::exists(spread_file)) {
      try {
        sub_disp_spread = DiskImageView<PixelMask<Vector2i> >(spread_file);
      }
      catch (...) {}
    }
  }

  ImageView<Matrix3x3> local_hom;
  if ( stereo_settings().seed_mode > 0 && stereo_settings().use_local_homography ){
    string local_hom_file = opt.out_prefix + "-local_hom.txt";
    read_local_homographies(local_hom_file, local_hom);
  }

  stereo::CostFunctionType cost_mode = get_cost_mode_value();
  Vector2i kernel_size    = stereo_settings().corr_kernel;
  BBox2i   trans_crop_win = stereo_settings().trans_crop_win;
  int      corr_timeout   = stereo_settings().corr_timeout;
  double   seconds_per_op = 0.0;
  if (corr_timeout > 0)
    seconds_per_op = calc_seconds_per_op(cost_mode, left_disk_image, right_disk_image, kernel_size);

  // Set up the reference to the stereo disparity code
  // - Processing is limited to trans_crop_win for use with parallel_stereo.
  ImageViewRef<PixelMask<Vector2f> > fullres_disparity =
    crop(SeededCorrelatorView( left_disk_image, right_disk_image, Lmask, Rmask,
                               sub_disp, sub_disp_spread, local_hom, kernel_size, 
                               cost_mode, corr_timeout, seconds_per_op),  
         trans_crop_win);

  // With SGM, we must do the entire image chunk as one tile. Otherwise,
  // if it gets done in smaller tiles, there will be artifacts at tile boundaries.
  bool using_sgm = (stereo_settings().stereo_algorithm > vw::stereo::CORRELATION_WINDOW);
  if (using_sgm) {
    Vector2i image_size = bounding_box(fullres_disparity).size();
    int max_dim = std::max(image_size[0], image_size[1]);
    if (stereo_settings().corr_tile_size_ovr < max_dim)
      vw_throw(ArgumentErr()
               << "Error: SGM processing is not permitted with a tile size smaller than the image!\n"
               << "Value of --corr-tile-size is " << stereo_settings().corr_tile_size_ovr
               << " but disparity size is " << image_size << ".\n" 
               << "Increase --corr-tile-size so the entire image fits in one tile, or "
               << "use parallel_stereo. Not that making --corr-tile-size larger than 9000 or so may "
               << "cause GDAL to crash.\n\n");
  }
  
  switch(stereo_settings().pre_filter_mode){
  case 2:
    vw_out() << "\t--> Using LOG pre-processing filter with "
             << stereo_settings().slogW << " sigma blur.\n"; 
    break;
  case 1:
    vw_out() << "\t--> Using Subtracted Mean pre-processing filter with "
	           << stereo_settings().slogW << " sigma blur.\n";
    break;
  default:
    vw_out() << "\t--> Using NO pre-processing filter." << endl;
  }

  cartography::GeoReference left_georef;
  bool   has_left_georef = read_georeference(left_georef,  opt.out_prefix + "-L.tif");
  bool   has_nodata      = false;
  double nodata          = -32768.0;

  string d_file = opt.out_prefix + "-D.tif";
  vw_out() << "Writing: " << d_file << "\n";
  if (stereo_settings().stereo_algorithm > vw::stereo::CORRELATION_WINDOW) {
    // SGM performs subpixel correlation in this step, so write out floats.
    vw::cartography::block_write_gdal_image(d_file, fullres_disparity,
			        has_left_georef, left_georef,
			        has_nodata, nodata, opt,
			        TerminalProgressCallback("asp", "\t--> Correlation :") );
  } else {
    // Otherwise cast back to integer results to save on storage space.
    vw::cartography::block_write_gdal_image(d_file, 
              pixel_cast<PixelMask<Vector2i> >(fullres_disparity),
			        has_left_georef, left_georef,
			        has_nodata, nodata, opt,
			        TerminalProgressCallback("asp", "\t--> Correlation :") );
  }
// Ricardo Monteiro - overwrite the homogrpahies
if ( stereo_settings().seed_mode > 0 && stereo_settings().use_local_homography ){
    string local_hom_file = opt.out_prefix + "-local_hom.txt";
    write_local_homographies(local_hom_file, local_hom);
    cout << "[Writing homographies]" << endl; 
  }

  vw_out() << "\n[ " << current_posix_time_string() << " ] : CORRELATION FINISHED \n";

} // End function stereo_correlation

int main(int argc, char* argv[]) {

  //try {
    xercesc::XMLPlatformUtils::Initialize();

    stereo_register_sessions();

    bool verbose = false;
    vector<ASPGlobalOptions> opt_vec;
    string output_prefix;
    asp::parse_multiview(argc, argv, CorrelationDescription(),
			 verbose, output_prefix, opt_vec);
    ASPGlobalOptions opt = opt_vec[0];

    // Leave the number of parallel block threads equal to the default unless we
    //  are using SGM in which case only one block at a time should be processed.
    // - Processing multiple blocks is possible, but it is better to use a larger blocks
    //   with more threads applied to the single block.
    // - Thread handling is still a little confusing because opt.num_threads is ONLY used
    //   to control the number of parallel image blocks written at a time.  Everything else
    //   reads directly from vw_settings().default_num_threads()
    const bool using_sgm = (stereo_settings().stereo_algorithm > vw::stereo::CORRELATION_WINDOW);
    opt.num_threads = vw_settings().default_num_threads();
    if (using_sgm)
      opt.num_threads = 1;

    // Integer correlator requires large tiles
    //---------------------------------------------------------
    int ts = stereo_settings().corr_tile_size_ovr;
    
    // GDAL block write sizes must be a multiple to 16 so if the input value is
    //  not a multiple of 16 increase it until it is.
    const int TILE_MULTIPLE = 16;
    if (ts % TILE_MULTIPLE != 0)
      ts = ((ts / TILE_MULTIPLE) + 1) * TILE_MULTIPLE;
      
    opt.raster_tile_size = Vector2i(ts, ts);

    // Internal Processes
    //---------------------------------------------------------
    stereo_correlation( opt );
  
    xercesc::XMLPlatformUtils::Terminate();
  //} ASP_STANDARD_CATCHES;

  return 0;
}

////// Ricardo Monteiro
BBox2f piecewiseAlignment_affineepipolar(	ImageView<float> left_image, 
										ImageView<float> right_image,
										ImageView<float> tile_left_image,
										ImageView<float> tile_right_image,	
										BBox2i bbox,
										Vector2i& left_size,
										Vector2i& right_size,
										vw::Matrix<double>& left_matrix,
										vw::Matrix<double>& right_matrix,
										BBox2f local_search_range)
{
	using namespace vw;

	double threshPiecewiseAlignment = 3;
	double avgDeltaY = -1.0;
	double threshRANSAC = 20.0;
	double threshSearchRange = 2;

	Matrix<double> H;
	double left_nodata_value  = numeric_limits<double>::quiet_NaN();
  	double right_nodata_value = numeric_limits<double>::quiet_NaN();
    std::vector<ip::InterestPoint> matched_ip1, matched_ip2;
	std::vector<ip::InterestPoint> matchedRANSAC_ip1,matchedRANSAC_ip2;
	std::vector<ip::InterestPoint> matchedRANSAC_final_ip1, matchedRANSAC_final_ip2;
	ip::InterestPoint aux_r_ip, aux_l_ip;
	bool success = false;
	char outputName[30]; // DEBUG
	int X = bbox.min().x()/ASPGlobalOptions::corr_tile_size();
	int Y = bbox.min().y()/ASPGlobalOptions::corr_tile_size();
	// detect and match ips
	sprintf(outputName, "matches_%d_%d", Y, X); // DEBUG
    try { success = homography_ip_matching1( tile_left_image, tile_right_image,
                                          stereo_settings().ip_per_tile,
                                          outputName, threshRANSAC, // before it was inlier_threshold
                                          left_nodata_value, right_nodata_value,
					  						matchedRANSAC_ip1, matchedRANSAC_ip2); }catch(...){}
	avgDeltaY = calcAverageDeltaY(matchedRANSAC_ip1, matchedRANSAC_ip2); // estimate global alignment
	cout << "[tile(" << Y << "," << X << ") avgDeltaY after global alignment = " << avgDeltaY << "]" << endl; // DEBUG
	if(avgDeltaY != -1 || avgDeltaY >= threshPiecewiseAlignment){ // if the alignment can be improved
		for ( size_t i = 0; i < matchedRANSAC_ip1.size(); i++ ) { // adjust ip matches to tile
				//cout << "[tile(" << Y << "," << X << " matchedRANSAC " << matchedRANSAC_ip2[i].y << " " << matchedRANSAC_ip2[i].x << "]" << endl; // DEBUG
	    	    //matchedRANSAC_ip1[i].x += bbox.min().x();
	    	    //matchedRANSAC_ip1[i].y += bbox.min().y();
	    	    //matchedRANSAC_ip2[i].x += bbox.min().x();
	    	    //matchedRANSAC_ip2[i].y += bbox.min().y();
		}
		sprintf(outputName, "matches_adj_%d_%d", Y, X); // DEBUG
		ip::write_binary_match_file(outputName, matchedRANSAC_ip1, matchedRANSAC_ip2); // DEBUG
		cout << "[tile(" << Y << "," << X << ")" << matchedRANSAC_ip1.size() << " matching points]" << endl;
		std::vector<Vector3> ransac_ip1 = iplist_to_vectorlist(matchedRANSAC_ip1), ransac_ip2 = iplist_to_vectorlist(matchedRANSAC_ip2);
		// RANSAC
		try {	
			left_size = affine_epipolar_rectification1(	left_size, right_size, matchedRANSAC_ip1, 												matchedRANSAC_ip2, left_matrix, right_matrix );
		} catch ( ... ) {
		  	left_matrix = math::identity_matrix<3>();
			right_matrix = math::identity_matrix<3>();
			return local_search_range;
		}
		// check left_matrix and right_matrix
		if(!check_homography_matrix(left_matrix, right_matrix, ransac_ip1, ransac_ip2, avgDeltaY, bbox)){
			left_matrix = math::identity_matrix<3>();
			right_matrix = math::identity_matrix<3>();
			return local_search_range;
		}
	}else{ // if the alignment cannot be improved
		left_matrix = math::identity_matrix<3>();
		right_matrix = math::identity_matrix<3>();
		return local_search_range;
	}
	return calcSearchRange(matchedRANSAC_ip1, matchedRANSAC_ip2, left_matrix, right_matrix, threshSearchRange);
}

vw::Matrix<double> piecewiseAlignment_homography(ImageView<float> left_image, 
			ImageView<float> right_image,
			ImageView<float> tile_left_image,
			ImageView<float> tile_right_image,	
			BBox2i bbox)
{
	using namespace vw;

	double threshPiecewiseAlignment = 3.0;
	double avgDeltaY = -1.0;
	double threshRANSAC = 3.0;

	Matrix<double> H;
	double left_nodata_value  = numeric_limits<double>::quiet_NaN();
  	double right_nodata_value = numeric_limits<double>::quiet_NaN();
    std::vector<ip::InterestPoint> matched_ip1, matched_ip2;
	std::vector<ip::InterestPoint> matchedRANSAC_ip1,matchedRANSAC_ip2;
	std::vector<ip::InterestPoint> matchedRANSAC_final_ip1, matchedRANSAC_final_ip2;
	ip::InterestPoint aux_r_ip, aux_l_ip;
	bool success = false;
	char outputName[30]; // DEBUG
	int X = bbox.min().x()/ASPGlobalOptions::corr_tile_size();
	int Y = bbox.min().y()/ASPGlobalOptions::corr_tile_size();
	// detect and match ips
	sprintf(outputName, "matches_%d_%d", Y, X); // DEBUG
    try { success = homography_ip_matching1( tile_left_image, tile_right_image,
                                          stereo_settings().ip_per_tile,
                                          outputName, threshRANSAC, // before it was inlier_threshold
                                          left_nodata_value, right_nodata_value,
					  						matchedRANSAC_ip1, matchedRANSAC_ip2); }catch(...){}

	avgDeltaY = calcAverageDeltaY(matchedRANSAC_ip1, matchedRANSAC_ip2); // estimate global alignment
	cout << "[tile(" << Y << "," << X << ") avgDeltaY after global alignment = " << avgDeltaY << "]" << endl; // DEBUG
	if(avgDeltaY != -1 || avgDeltaY >= threshPiecewiseAlignment){ // if the alignment can be improved
		for ( size_t i = 0; i < matchedRANSAC_ip1.size(); i++ ) { // adjust ip matches to tile
				//cout << "[tile(" << Y << "," << X << " matchedRANSAC " << matchedRANSAC_ip2[i].y << " " << matchedRANSAC_ip2[i].x << "]" << endl; // DEBUG
	    	    matchedRANSAC_ip1[i].x += bbox.min().x();
	    	    matchedRANSAC_ip1[i].y += bbox.min().y();
	    	    matchedRANSAC_ip2[i].x += bbox.min().x();
	    	    matchedRANSAC_ip2[i].y += bbox.min().y();
		}
		sprintf(outputName, "matches_adj_%d_%d", Y, X); // DEBUG
		ip::write_binary_match_file(outputName, matchedRANSAC_ip1, matchedRANSAC_ip2); // DEBUG
		cout << "[tile(" << Y << "," << X << ")" << matchedRANSAC_ip1.size() << " matching points]" << endl;
		std::vector<Vector3> ransac_ip1 = iplist_to_vectorlist(matchedRANSAC_ip1), ransac_ip2 = iplist_to_vectorlist(matchedRANSAC_ip2);
		std::vector<size_t> indices;
		// RANSAC
		try {
			// use RANSAC to try to fit a homography with the ip matches
			typedef math::RandomSampleConsensus<math::HomographyFittingFunctor, math::InterestPointErrorMetric> RansacT;
		  	const int    MIN_NUM_OUTPUT_INLIERS = ransac_ip1.size()/2;
		   // const int    MIN_NUM_OUTPUT_INLIERS = 4;
		  	const int    NUM_ITERATIONS         = 200;
		  	RansacT ransac( math::HomographyFittingFunctor(), math::InterestPointErrorMetric(), NUM_ITERATIONS,
				  threshRANSAC, MIN_NUM_OUTPUT_INLIERS, true );
		  	H = ransac(ransac_ip2,ransac_ip1); // 2 then 1 is used here for legacy reasons
			indices = ransac.inlier_indices(H,ransac_ip2,ransac_ip1);
			BOOST_FOREACH( size_t& index, indices ){
				aux_l_ip.x = ransac_ip1[index].x(); // DEBUG
    			aux_l_ip.y = ransac_ip1[index].y(); // DEBUG
    			aux_r_ip.x = ransac_ip2[index].x(); // DEBUG
    			aux_r_ip.y = ransac_ip2[index].y(); // DEBUG
				matchedRANSAC_final_ip1.push_back(aux_r_ip); // DEBUG
				matchedRANSAC_final_ip2.push_back(aux_l_ip); // DEBUG
			}
			cout << "[tile(" << Y << "," << X << ")" << matchedRANSAC_final_ip1.size() << " matching points after H]" << endl;	
		} catch ( ... ) {
		  	return math::identity_matrix<3>();
		}

		// check H
		//if(check_homography_matrix(H, iplist_to_vectorlist(matchedRANSAC_final_ip1), iplist_to_vectorlist(matchedRANSAC_final_ip2), indices, avgDeltaY, bbox))
		if(check_homography_matrix(H, ransac_ip1, ransac_ip2, indices, avgDeltaY, bbox))
			//return H; // good H
			return math::HomographyFittingFunctor()(ransac_ip2, ransac_ip1, H);
		else
			return math::identity_matrix<3>(); // bad H
	}else // if the alignment cannot be improved
		return math::identity_matrix<3>();
}



vw::Matrix<double> piecewiseAlignment(ImageView<float> left_image, 
			ImageView<float> right_image,
			ImageView<float> tile_left_image,
			ImageView<float> tile_right_image,	
			BBox2i bbox)
{
	vw::Matrix<double> fullres_hom;
	float left_nodata_value  = numeric_limits<float>::quiet_NaN();
  	float right_nodata_value = numeric_limits<float>::quiet_NaN();
	std::vector<ip::InterestPoint> left_ip, right_ip;
	bool success = false;
	char outputName[30];
	int ts = ASPGlobalOptions::corr_tile_size();
	int W = bbox.min().x()/ts;
	int H = bbox.min().y()/ts;
	double avgDeltaY = -1.0/*, avgDeltaY_afterPiecewise = -1.0*/;
	double threshRANSAC = 1.0;
	double threshPiecewiseAlignment = 3.0;
	sprintf(outputName, "matches_%d_%d", H, W);
	try { success = homography_ip_matching1( tile_left_image, tile_right_image,
                                          stereo_settings().ip_per_tile,
                                          outputName, threshRANSAC, // before it was inlier_threshold
                                          left_nodata_value, right_nodata_value,
					  left_ip, right_ip); }catch(...){}
	cout << "[tile(" << H << "," << W << ")" << left_ip.size() << " matching points]" << endl;
        //for(size_t i = 0; i < right_ip.size(); i++)
           //cout << "ip matches after ip match with RANSAC " << right_ip[i].y << "\n";
	avgDeltaY = calcAverageDeltaY(left_ip, right_ip);
	if(avgDeltaY == -1 || avgDeltaY < threshPiecewiseAlignment)
	   success = false;
	cout << "[tile(" << H << "," << W << ") avgDeltaY after global alignment = " << avgDeltaY << "]" << endl;
	if(success)
	{
	    cout << "[tile(" << H << "," << W << ") success!]" << endl;
	    Matrix<double> left_matrix = math::identity_matrix<3>();
	    Matrix<double> right_matrix = fullres_hom;
	    try {
		// adjust the ip to the resolution of the full image
		for ( size_t i = 0; i < left_ip.size(); i++ ) {
	    	    left_ip[i].x += bbox.min().x();
	    	    left_ip[i].y += bbox.min().y();
	    	    right_ip[i].x += bbox.min().x();
	    	    right_ip[i].y += bbox.min().y();
		}
		sprintf(outputName, "matches_adj_%d_%d.tif", H, W);
		 
		ip::write_binary_match_file(outputName, left_ip, right_ip); // write ip matches after adjustment
	        homography_rectification1( false,
                                left_image.get_size(), right_image.get_size(),
                                left_ip, right_ip, left_matrix, right_matrix, 
				threshRANSAC, avgDeltaY, bbox );
		// // if true use:
		//right_matrix(0,2) -= left_matrix(0,2);
      		//right_matrix(1,2) -= left_matrix(1,2); 
		fullres_hom = right_matrix; // overwrite fullres_hom
		cout << "[tile(" << H << "," << W << ") updated fullres_hom]" << endl; 
	    }
    	    catch ( ... ){
		fullres_hom = math::identity_matrix<3>(); // overwrite fullres_hom
		cout << "[tile(" << H << "," << W << ") updated fullres_hom with identity matrix]" << endl; 
	    }
	    cout << "[tile(" << H << "," << W << ") " << fullres_hom << "]" << endl;
	}else{
	    cout << "[tile(" << H << "," << W << ") NO success!]" << endl;
	    fullres_hom = math::identity_matrix<3>(); // overwrite fullres_hom
	    cout << "[tile(" << H << "," << W << ") updated fullres_hom with identity matrix]" << endl; 
	}

	return fullres_hom;
}

Vector2i
  affine_epipolar_rectification1( Vector2i const& left_size,
                                 Vector2i const& right_size,
                                 std::vector<ip::InterestPoint> const& ip1,
                                 std::vector<ip::InterestPoint> const& ip2,
                                 Matrix<double>& left_matrix,
                                 Matrix<double>& right_matrix ) {
    // Create affine fundamental matrix
    Matrix<double> fund = linear_affine_fundamental_matrix( ip1, ip2 );

    // Solve for rotation matrices
    double Hl = sqrt( fund(2,0)*fund(2,0) + fund(2,1)*fund(2,1) );
    double Hr = sqrt( fund(0,2)*fund(0,2) + fund(1,2)*fund(1,2) );
    Vector2 epipole(-fund(2,1),fund(2,0)), epipole_prime(-fund(1,2),fund(0,2));
    if ( epipole.x() < 0 )
      epipole = -epipole;
    if ( epipole_prime.x() < 0 )
      epipole_prime = -epipole_prime;
    epipole.y() = -epipole.y();
    epipole_prime.y() = -epipole_prime.y();

    left_matrix = math::identity_matrix<3>();
    right_matrix = math::identity_matrix<3>();
    left_matrix(0,0) = epipole[0]/Hl;
    left_matrix(0,1) = -epipole[1]/Hl;
    left_matrix(1,0) = epipole[1]/Hl;
    left_matrix(1,1) = epipole[0]/Hl;
    right_matrix(0,0) = epipole_prime[0]/Hr;
    right_matrix(0,1) = -epipole_prime[1]/Hr;
    right_matrix(1,0) = epipole_prime[1]/Hr;
    right_matrix(1,1) = epipole_prime[0]/Hr;

    // Solve for ideal scaling and translation
    solve_y_scaling( ip1, ip2, left_matrix, right_matrix );

    // Solve for ideal shear, scale, and translation of X axis
    solve_x_shear( ip1, ip2, left_matrix, right_matrix );

    // Work out the ideal render size.
    BBox2i output_bbox, right_bbox;
    output_bbox.grow( subvector(left_matrix*Vector3(0,0,1),0,2) );
    output_bbox.grow( subvector(left_matrix*Vector3(left_size.x(),0,1),0,2) );
    output_bbox.grow( subvector(left_matrix*Vector3(left_size.x(),left_size.y(),1),0,2) );
    output_bbox.grow( subvector(left_matrix*Vector3(0,left_size.y(),1),0,2) );
    right_bbox.grow( subvector(right_matrix*Vector3(0,0,1),0,2) );
    right_bbox.grow( subvector(right_matrix*Vector3(right_size.x(),0,1),0,2) );
    right_bbox.grow( subvector(right_matrix*Vector3(right_size.x(),right_size.y(),1),0,2) );
    right_bbox.grow( subvector(right_matrix*Vector3(0,right_size.y(),1),0,2) );
    output_bbox.crop( right_bbox );

    left_matrix(0,2) -= output_bbox.min().x();
    right_matrix(0,2) -= output_bbox.min().x();
    left_matrix(1,2) -= output_bbox.min().y();
    right_matrix(1,2) -= output_bbox.min().y();

    return Vector2i( output_bbox.width(), output_bbox.height() );
  }



  // Homography IP matching - Ricardo Monteiro - return ip matching
  //
  // This applies only the homography constraint. Not the best...
  template <class Image1T, class Image2T>
  bool homography_ip_matching1( vw::ImageViewBase<Image1T> const& image1,
			       vw::ImageViewBase<Image2T> const& image2,
			       int ip_per_tile,
			       std::string const& output_name,
			       int inlier_threshold,
			       double nodata1,
			       double nodata2,
			       std::vector<ip::InterestPoint>& final_ip1,
 			       std::vector<ip::InterestPoint>& final_ip2) {

    using namespace vw;

    std::vector<ip::InterestPoint> matched_ip1, matched_ip2;
    detect_match_ip( matched_ip1, matched_ip2,
		     image1.impl(), image2.impl(),
		     ip_per_tile,
		     nodata1, nodata2 );
    cout << "matches left = " <<  matched_ip1.size() << " matches right = " <<  matched_ip2.size() << endl;
    if ( matched_ip1.size() == 0 || matched_ip2.size() == 0 )
      return false;
    std::vector<Vector3> ransac_ip1 = iplist_to_vectorlist(matched_ip1),
			 			 ransac_ip2 = iplist_to_vectorlist(matched_ip2);
    std::vector<size_t> indices;
    try {
      typedef math::RandomSampleConsensus<math::HomographyFittingFunctor, math::InterestPointErrorMetric> RansacT;
      const int    MIN_NUM_OUTPUT_INLIERS = ransac_ip1.size()/2;
     //const int    MIN_NUM_OUTPUT_INLIERS = ransac_ip1.size()/10;
	 // const int    MIN_NUM_OUTPUT_INLIERS = 4;
      const int    NUM_ITERATIONS         = 100;
      RansacT ransac( math::HomographyFittingFunctor(),
		      math::InterestPointErrorMetric(), NUM_ITERATIONS,
		      inlier_threshold,
		      MIN_NUM_OUTPUT_INLIERS, true
		      );
      
      Matrix<double> H(ransac(ransac_ip2,ransac_ip1)); // 2 then 1 is used here for legacy reasons
      //vw_out() << "\t--> Homography: " << H << "\n";
     // cout << "homography_ip_matching " << H << " ";
      indices = ransac.inlier_indices(H,ransac_ip2,ransac_ip1);
     // cout << H << endl;
    } catch (const math::RANSACErr& e ) {
      //vw_out() << "RANSAC Failed: " << e.what() << "\n";
      return false;
    }

   // std::vector<ip::InterestPoint> final_ip1, final_ip2;
int i = 0; // DEBUG
    BOOST_FOREACH( size_t& index, indices ) {
      final_ip1.push_back(matched_ip1[index]);
      final_ip2.push_back(matched_ip2[index]);
//cout << " after ip matching  " << final_ip2[i].y << " " << final_ip2[i].x << "]" << endl; // DEBUG
i++;  // DEBUG
    }


    //// DEBUG - Draw out the point matches pre-geometric filtering
    //vw_out() << "\t    Writing IP debug image2! " << std::endl;
    //write_match_image("InterestPointMatching__ip_matching_debug2.tif",
    //                  image1, image2,
    //                  final_ip1, final_ip2);

    //vw_out() << "\t    * Writing match file: " << output_name << "\n";
    ip::write_binary_match_file(output_name, final_ip1, final_ip2);
    return true;
  }

Vector2i
  homography_rectification1( bool adjust_left_image_size,
			    Vector2i const& left_size,
			    Vector2i const& right_size,
			    std::vector<ip::InterestPoint> const& left_ip,
			    std::vector<ip::InterestPoint> const& right_ip,
			    vw::Matrix<double>& left_matrix,
			    vw::Matrix<double>& right_matrix,
			    double threshRANSAC,
			    double minAvgDeltaY, BBox2i bbox ) {
    // Reformat the interest points for RANSAC
    std::vector<Vector3>  right_copy = iplist_to_vectorlist(right_ip),
			  left_copy  = iplist_to_vectorlist(left_ip);

    //double thresh_factor = stereo_settings().ip_inlier_factor; // 1/15 by default
    
    // Use RANSAC to determine a good homography transform between the images
    math::RandomSampleConsensus<math::HomographyFittingFunctor, math::InterestPointErrorMetric>
      ransac( math::HomographyFittingFunctor(),
	      math::InterestPointErrorMetric(),
	      100, // num iter
              //threshRANSAC * norm_2(Vector2(left_size.x(),left_size.y())) * (1.5*thresh_factor), // inlier thresh 
	      threshRANSAC, // Ricardo Monteiro //////////
	      left_copy.size()*1/10 // min output inliers
	      );

  //  std::cout << "[RANSAC old thresh = " << norm_2(Vector2(left_size.x(),left_size.y())) * (1.5*thresh_factor) << "]\n";
  //  std::cout << "[RANSAC new thresh = " << threshRANSAC << "]\n";
	
    Matrix<double> H = ransac(right_copy, left_copy);
    cout << "homography_rectification " << H << " ";
    std::vector<size_t> indices = ransac.inlier_indices(H, right_copy, left_copy);
    cout << H;

    if(check_homography_matrix(H, left_copy, right_copy, indices, minAvgDeltaY, bbox)){
    // Set right to a homography that has been refined just to our inliers
    	left_matrix  = math::identity_matrix<3>();
    	right_matrix = math::HomographyFittingFunctor()(right_copy, left_copy, H);
    cout << H << endl;;
    }else{
	left_matrix  = math::identity_matrix<3>();
    	right_matrix = math::identity_matrix<3>();
    }

    // Work out the ideal render size
    BBox2i output_bbox, right_bbox;
    output_bbox.grow( Vector2i(0,0) );
    output_bbox.grow( Vector2i(left_size.x(),0) );
    output_bbox.grow( Vector2i(0,left_size.y()) );
    output_bbox.grow( left_size );

    if (adjust_left_image_size){
      // Crop the left and right images to the shared region. This is
      // done for efficiency.  It may not be always desirable though,
      // as in this case we lose the one-to-one correspondence between
      // original input left image pixels and output disparity/point
      // cloud pixels.
      Vector3 temp = right_matrix*Vector3(0,0,1);
      temp /= temp.z();
      right_bbox.grow( subvector(temp,0,2) );
      temp = right_matrix*Vector3(right_size.x(),0,1);
      temp /= temp.z();
      right_bbox.grow( subvector(temp,0,2) );
      temp = right_matrix*Vector3(0,right_size.y(),1);
      temp /= temp.z();
      right_bbox.grow( subvector(temp,0,2) );
      temp = right_matrix*Vector3(right_size.x(),right_size.y(),1);
      temp /= temp.z();
      right_bbox.grow( subvector(temp,0,2) );

      output_bbox.crop( right_bbox );

      //  Move the ideal render size to be aligned up with origin
      left_matrix (0,2) -= output_bbox.min().x();
      right_matrix(0,2) -= output_bbox.min().x();
      left_matrix (1,2) -= output_bbox.min().y();
      right_matrix(1,2) -= output_bbox.min().y();
    }

    return Vector2i( output_bbox.width(), output_bbox.height() );
  }


bool check_homography_matrix(	Matrix<double>       const& left_matrix,
								Matrix<double>       const& right_matrix,
			       				std::vector<Vector3> const& left_points,
			       				std::vector<Vector3> const& right_points,
			       				double minAvgDeltaY, 
								BBox2i bbox
								){

    // Sanity checks. If these fail, most likely the two images are too different
    // for stereo to succeed.
    /*if ( indices.size() < std::min( right_points.size(), left_points.size() )/2 ){
      vw_out(WarningMessage) << "InterestPointMatching: The number of inliers is less "
                             << "than 1/2 of the number of points. The inputs may be invalid.\n";
	return false;
    }*/

    double det = fabs(left_matrix(0, 0)*left_matrix(1, 1) - left_matrix(0, 1)*left_matrix(1, 0));
    if (det <= 0.1 || det >= 10.0){
      vw_out(WarningMessage) << "InterestPointMatching: The determinant of the 2x2 submatrix "
                             << "of the homography matrix " << left_matrix << " is " << det
                             << ". There could be a large scale discrepancy among the input images "
                             << "or the inputs may be an invalid stereo pair.\n";
	return false;
    }
	det = fabs(right_matrix(0, 0)*right_matrix(1, 1) - right_matrix(0, 1)*right_matrix(1, 0));
    if (det <= 0.1 || det >= 10.0){
      vw_out(WarningMessage) << "InterestPointMatching: The determinant of the 2x2 submatrix "
                             << "of the homography matrix " << right_matrix << " is " << det
                             << ". There could be a large scale discrepancy among the input images "
                             << "or the inputs may be an invalid stereo pair.\n";
	return false;
    }

    // check if the avgDeltaY after piecewise alignment is better than the minAvgDeltaY
    std::vector<Vector3> right_ip; 
    std::vector<Vector3> left_ip;
	std::vector<ip::InterestPoint> r_ip, l_ip;
	ip::InterestPoint aux_r_ip, aux_l_ip;
    double avgDeltaY = -1;
    int ts = ASPGlobalOptions::corr_tile_size();
    for(size_t i = 0; i < right_points.size(); i++)
    { 
		//cout << " ip matchings " << right_points[i].y() << " " << right_points[i].x() << "]" << endl; // DEBUG
        right_ip.push_back(right_matrix * Vector3(right_points[i].x(), right_points[i].y(), 1));
        left_ip.push_back(left_matrix * Vector3(left_points[i].x(), left_points[i].y(), 1));
	// Normalize the coordinates, but don't divide by 0
        if (right_ip[i].z() == 0 || left_ip[i].z() == 0) 
            continue;
        right_ip[i] /= right_ip[i].z();
        left_ip[i] /= left_ip[i].z();

		aux_l_ip.x = left_ip[i].x() /*- bbox.min().x()*/; // DEBUG
    	aux_l_ip.y = left_ip[i].y() /*- bbox.min().y()*/; // DEBUG
    	aux_r_ip.x = right_ip[i].x() /*- bbox.min().x()*/; // DEBUG
    	aux_r_ip.y = right_ip[i].y() /*- bbox.min().y()*/; // DEBUG
		r_ip.push_back(aux_r_ip); // DEBUG
		l_ip.push_back(aux_l_ip); // DEBUG
		//cout << " ip matchings after H " << right_ip[i].y() << " " << right_ip[i].x() << "]" << endl; // DEBUG
    }
    avgDeltaY = calcAverageDeltaY(left_ip, right_ip);
    cout << "[tile(" << bbox.min().y()/ts << "," << bbox.min().x()/ts << ") avgDeltaY after piecewise alignment = " << avgDeltaY << "]" << endl;
    
	char outputName[30]; // DEBUG
	int X = bbox.min().x()/ASPGlobalOptions::corr_tile_size(); // DEBUG
	int Y = bbox.min().y()/ASPGlobalOptions::corr_tile_size(); // DEBUG
	sprintf(outputName, "matches_after_H_%d_%d", Y, X); // DEBUG
	ip::write_binary_match_file(outputName, l_ip, r_ip); // DEBUG
    if(avgDeltaY == -1 || avgDeltaY >= minAvgDeltaY)
        return false;

    return true;

  }

bool check_homography_matrix(Matrix<double>       const& H,
			       std::vector<Vector3> const& left_points,
			       std::vector<Vector3> const& right_points,
			       std::vector<size_t>  const& indices,
			       double minAvgDeltaY, BBox2i bbox
			       ){

    // Sanity checks. If these fail, most likely the two images are too different
    // for stereo to succeed.
    /*if ( indices.size() < std::min( right_points.size(), left_points.size() )/2 ){
      vw_out(WarningMessage) << "InterestPointMatching: The number of inliers is less "
                             << "than 1/2 of the number of points. The inputs may be invalid.\n";
	return false;
    }*/

    double det = fabs(H(0, 0)*H(1, 1) - H(0, 1)*H(1, 0));
    if (det <= 0.1 || det >= 10.0){
      vw_out(WarningMessage) << "InterestPointMatching: The determinant of the 2x2 submatrix "
                             << "of the homography matrix " << H << " is " << det
                             << ". There could be a large scale discrepancy among the input images "
                             << "or the inputs may be an invalid stereo pair.\n";
	return false;
    }

    // check if the avgDeltaY after piecewise alignment is better than the minAvgDeltaY
    std::vector<Vector3> right_ip; 
    std::vector<Vector3> left_ip;
	std::vector<ip::InterestPoint> r_ip, l_ip;
	ip::InterestPoint aux_r_ip, aux_l_ip;
    double avgDeltaY = -1;
    int ts = ASPGlobalOptions::corr_tile_size();
    for(size_t i = 0; i < right_points.size(); i++)
    { 
		//cout << " ip matchings " << right_points[i].y() << " " << right_points[i].x() << "]" << endl; // DEBUG
        right_ip.push_back(H * Vector3(right_points[i].x(), right_points[i].y(), 1));
        left_ip.push_back(Vector3(left_points[i].x(), left_points[i].y(), 1));
	// Normalize the coordinates, but don't divide by 0
        if (right_ip[i].z() == 0 || left_ip[i].z() == 0) 
            continue;
        right_ip[i] /= right_ip[i].z();
        left_ip[i] /= left_ip[i].z();

		aux_l_ip.x = left_ip[i].x() - bbox.min().x(); // DEBUG
    	aux_l_ip.y = left_ip[i].y() - bbox.min().y(); // DEBUG
    	aux_r_ip.x = right_ip[i].x() - bbox.min().x(); // DEBUG
    	aux_r_ip.y = right_ip[i].y() - bbox.min().y(); // DEBUG
		r_ip.push_back(aux_r_ip); // DEBUG
		l_ip.push_back(aux_l_ip); // DEBUG
		//cout << " ip matchings after H " << right_ip[i].y() << " " << right_ip[i].x() << "]" << endl; // DEBUG
    }
    avgDeltaY = calcAverageDeltaY(left_ip, right_ip);
    cout << "[tile(" << bbox.min().y()/ts << "," << bbox.min().x()/ts << ") avgDeltaY after piecewise alignment = " << avgDeltaY << "]" << endl;
    
	char outputName[30]; // DEBUG
	int X = bbox.min().x()/ASPGlobalOptions::corr_tile_size(); // DEBUG
	int Y = bbox.min().y()/ASPGlobalOptions::corr_tile_size(); // DEBUG
	sprintf(outputName, "matches_after_H_%d_%d", Y, X); // DEBUG
	ip::write_binary_match_file(outputName, l_ip, r_ip); // DEBUG
    if(avgDeltaY == -1 || avgDeltaY >= minAvgDeltaY)
        return false;

    return true;

  }

double calcAverageDeltaY(std::vector<ip::InterestPoint> const& left_points, std::vector<ip::InterestPoint> const& right_points)
{
    double accuDiff = 0;

    if(left_points.size()){
        for ( size_t i = 0; i < left_points.size(); i++ )
	    accuDiff += abs(left_points[i].y - right_points[i].y);
        return accuDiff/left_points.size(); // average
    }else
	return -1; // not valid
}

double calcAverageDeltaY(std::vector<Vector3> const& left_points, std::vector<Vector3> const& right_points)
{
    double accuDiff = 0;

    if(left_points.size()){
        for ( size_t i = 0; i < left_points.size(); i++ )
	    accuDiff += abs(left_points[i].y() - right_points[i].y());
        return accuDiff/left_points.size(); // average
    }else
	return -1; // not valid
}

BBox2f calcSearchRange(std::vector<ip::InterestPoint> const& left_ip, std::vector<ip::InterestPoint> const& right_ip, Matrix<double> const& left_matrix, Matrix<double> const& right_matrix, double multi)
{
	std::vector<int> diffY, diffX;
	int maxDiffY = 0, maxDiffX = 0, minDiffY = 0, minDiffX = 0; 
	std::vector<Vector3> trans_left_points, trans_right_points;

	for(size_t i = 0; i < left_ip.size(); i++) // transform ip matches
    { 
        trans_right_points.push_back(right_matrix * Vector3(right_ip[i].x, right_ip[i].y, 1));
        trans_left_points.push_back(left_matrix * Vector3(left_ip[i].x, left_ip[i].y, 1));
        if (trans_right_points[i].z() == 0 || trans_left_points[i].z() == 0) 
            continue;
        trans_right_points[i] /= trans_right_points[i].z();
        trans_left_points[i] /= trans_left_points[i].z();
    }

	for ( size_t i = 0; i < trans_right_points.size(); i++ ){ // gen list of diff (right - left)
		diffY.push_back(trans_right_points[i].y() - trans_left_points[i].y());
		diffX.push_back(trans_right_points[i].x() - trans_left_points[i].x());
	}
	for ( size_t i = 0; i < trans_right_points.size(); i++ ){ // get min and max
		if(diffY[i] < minDiffY)
			minDiffY = diffY[i];
		if(diffX[i] < minDiffX)
			minDiffX = diffX[i];
		if(diffY[i] > maxDiffY)
			maxDiffY = diffY[i];
		if(diffX[i] > maxDiffX)
			maxDiffX = diffX[i];
	}
	return BBox2f(multi * minDiffX, multi * minDiffY, (multi * maxDiffX) - (multi * minDiffX), (multi * maxDiffY) - (multi * minDiffY));
}

