# Set up environments for
# 1. fsl (tested with 6.0.1)
# 2. freesurfer (tested with 5.3.0)
# 3. matlab (tested with R2018a)
source /cluster/software/fsl/6.0.1/setup.source
module load freesurfer/5.3.0
export FSLOUTPUTTYPE="NIFTI_GZ"
module load matlab/R2018a

# Declare 
# 1. the name of the subject to analyse: s
# 2. the directory that contains the subject: subdir
# 3. the directory to which the output is stored: outdir
# 4. the directory that contains the matlab scripts: scriptdir
s=bert
subdir=/path/to/data
outdir=/path/to/results
scriptdir=/path/to/scripts

# initialize 
export SUBJECTS_DIR=${subdir}
hemi="lh"
cd $outdir

#run watershed
mri_watershed -surf bemsurf -brainsurf brainsurf -useSRAS ${subdir}/${s}/mri/nu.mgz watershedbrain.mgz
cp lh.brainsurf* ${subdir}/${s}/surf/

#convert the nu.mgz to nifti                                                                                                                                                                                              
mri_convert ${subdir}/${s}/mri/nu.mgz ${outdir}/nu.nii.gz
mri_convert ${subdir}/${s}/mri/aseg.mgz ${outdir}/aseg.nii.gz

# make individual mask
matlab -nodesktop -nosplash -r "addpath ${scriptdir}; makeIndividualMask $s; exit"
mri_convert ${outdir}/${s}_individual_mask.mgz ${outdir}/${s}_individual_mask.nii.gz
mask="${outdir}/${s}_individual_mask.nii.gz"

#now multiply with the mask
fslmaths ${outdir}/nu.nii.gz -mul ${mask} ${outdir}/nu_masked.nii.gz

# and reorient
fslreorient2std ${outdir}/nu.nii.gz ${outdir}/nu_reorient.nii.gz
fslreorient2std ${outdir}/nu_masked.nii.gz ${outdir}/nu_masked_reorient.nii.gz

#convert that back to mgz format
mri_convert  ${outdir}/nu_masked_reorient.nii.gz  ${outdir}/nu_masked_reorient.mgz

#run stepwise analysis
for step in $(seq 0.5 0.5 25)
 do
  mri_vol2surf \
  --mov ${outdir}/nu_masked_reorient.mgz \
  --hemi $hemi --noreshape --interp trilinear \
  --surf brainsurf_outer_skin_surface \
  --projdist -${step} \
  --srchit ${outdir}/steps_${step}mm.mgz \
  --o ${outdir}/steps_${step}mm.nu.mgh --regheader $s

  mri_convert ${outdir}/steps_${step}mm.mgz ${outdir}/steps_${step}mm.nii.gz
done

#run matlab analysis
matlab -nodesktop -nosplash -r "addpath ${scriptdir}; compileLayers $s; exit"

# now there is a *_layers.csv file that can be fed into the neural net model
