%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 3944 $)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.temporal.st.scans = {'<UNDEFINED>'};
matlabbatch{1}.spm.temporal.st.nslices = 46;
matlabbatch{1}.spm.temporal.st.tr = 3;
matlabbatch{1}.spm.temporal.st.ta = 2.93478260869565;
matlabbatch{1}.spm.temporal.st.so = [2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45];
matlabbatch{1}.spm.temporal.st.refslice = 23;
matlabbatch{1}.spm.temporal.st.prefix = 'a';
matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep;
matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1).tname = 'Session';
matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1).tgt_spec = {};
matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1).sname = 'Slice Timing: Slice Timing Corr. Images (Sess 1)';
matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1).src_output = substruct('()',{1}, '.','files');
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = {''};
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [0 1];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
matlabbatch{3}.spm.spatial.coreg.estimate.ref(1) = cfg_dep;
matlabbatch{3}.spm.spatial.coreg.estimate.ref(1).tname = 'Reference Image';
matlabbatch{3}.spm.spatial.coreg.estimate.ref(1).tgt_spec = {};
matlabbatch{3}.spm.spatial.coreg.estimate.ref(1).sname = 'Realign: Estimate & Reslice: Mean Image';
matlabbatch{3}.spm.spatial.coreg.estimate.ref(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.spm.spatial.coreg.estimate.ref(1).src_output = substruct('.','rmean');
matlabbatch{3}.spm.spatial.coreg.estimate.source = '<UNDEFINED>';
matlabbatch{3}.spm.spatial.coreg.estimate.other = '<UNDEFINED>';
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
matlabbatch{4}.spm.spatial.preproc.data = '<UNDEFINED>';
matlabbatch{4}.spm.spatial.preproc.output.GM = [0 0 1];
matlabbatch{4}.spm.spatial.preproc.output.WM = [0 0 1];
matlabbatch{4}.spm.spatial.preproc.output.CSF = [0 0 1];
matlabbatch{4}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{4}.spm.spatial.preproc.output.cleanup = 0;
matlabbatch{4}.spm.spatial.preproc.opts.tpm = {
                                               'C:\Matlab_toolboxes\spm8\tpm\grey.nii'
                                               'C:\Matlab_toolboxes\spm8\tpm\white.nii'
                                               'C:\Matlab_toolboxes\spm8\tpm\csf.nii'
                                               };
matlabbatch{4}.spm.spatial.preproc.opts.ngaus = [2
                                                 2
                                                 2
                                                 4];
matlabbatch{4}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{4}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{4}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{4}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{4}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{4}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{4}.spm.spatial.preproc.opts.msk = '<UNDEFINED>';
matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1) = cfg_dep;
matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1).tname = 'Parameter File';
matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1).tgt_spec = {};
matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1).sname = 'Segment: Norm Params Subj->MNI';
matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{5}.spm.spatial.normalise.write.subj.matname(1).src_output = substruct('()',{1}, '.','snfile', '()',{':'});
matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep;
matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1).tname = 'Images to Write';
matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1).tgt_spec = {};
matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1).sname = 'Realign: Estimate & Reslice: Realigned Images (Sess 1)';
matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1).src_output = substruct('.','sess', '()',{1}, '.','cfiles');
matlabbatch{5}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{5}.spm.spatial.normalise.write.roptions.bb = [-91 -125 -71
                                                          89 91 109];
matlabbatch{5}.spm.spatial.normalise.write.roptions.vox = [3 3 3];
matlabbatch{5}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{5}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{5}.spm.spatial.normalise.write.roptions.prefix = 'w';
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1) = cfg_dep;
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1).tname = 'Parameter File';
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1).tgt_spec = {};
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1).sname = 'Segment: Norm Params Subj->MNI';
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1).src_output = substruct('()',{1}, '.','snfile', '()',{':'});
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep;
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1).tname = 'Images to Write';
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1).tgt_spec = {};
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1).sname = 'Segment: Bias Corr Images';
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1).src_output = substruct('()',{1}, '.','biascorr', '()',{':'});
matlabbatch{6}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{6}.spm.spatial.normalise.write.roptions.bb = [-91 -125 -71
                                                          89 91 109];
matlabbatch{6}.spm.spatial.normalise.write.roptions.vox = [1 1 1];
matlabbatch{6}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{6}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{6}.spm.spatial.normalise.write.roptions.prefix = 'w';
matlabbatch{7}.spm.spatial.normalise.write.subj.matname(1) = cfg_dep;
matlabbatch{7}.spm.spatial.normalise.write.subj.matname(1).tname = 'Parameter File';
matlabbatch{7}.spm.spatial.normalise.write.subj.matname(1).tgt_spec = {};
matlabbatch{7}.spm.spatial.normalise.write.subj.matname(1).sname = 'Segment: Norm Params Subj->MNI';
matlabbatch{7}.spm.spatial.normalise.write.subj.matname(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{7}.spm.spatial.normalise.write.subj.matname(1).src_output = substruct('()',{1}, '.','snfile', '()',{':'});
matlabbatch{7}.spm.spatial.normalise.write.subj.resample = '<UNDEFINED>';
matlabbatch{7}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{7}.spm.spatial.normalise.write.roptions.bb = [-91 -125 -71
                                                          89 91 109];
matlabbatch{7}.spm.spatial.normalise.write.roptions.vox = [1 1 1];
matlabbatch{7}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{7}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{7}.spm.spatial.normalise.write.roptions.prefix = 'w';
matlabbatch{8}.spm.spatial.smooth.data(1) = cfg_dep;
matlabbatch{8}.spm.spatial.smooth.data(1).tname = 'Images to Smooth';
matlabbatch{8}.spm.spatial.smooth.data(1).tgt_spec = {};
matlabbatch{8}.spm.spatial.smooth.data(1).sname = 'Normalise: Write: Normalised Images (Subj 1)';
matlabbatch{8}.spm.spatial.smooth.data(1).src_exbranch = substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{8}.spm.spatial.smooth.data(1).src_output = substruct('()',{1}, '.','files');
matlabbatch{8}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{8}.spm.spatial.smooth.dtype = 0;
matlabbatch{8}.spm.spatial.smooth.im = 0;
matlabbatch{8}.spm.spatial.smooth.prefix = 's';
