using Test
using SCTransform
using SCTransform: scparams_estimate, scparams_detect_outliers, scparams_bandwidth, scparams_regularize
using SingleCell10x

using DelimitedFiles
using DataFrames


function read_model_pars(fn, delim=',')
	table,header = readdlm(fn, delim; header=true)
	@assert header==["" "theta" "beta0" "beta1"]

	genename = String.(table[:,1])
	theta = Float64.(table[:,2])
	beta0 = Float64.(table[:,3])
	beta1 = Float64.(table[:,4])
	(;genename, theta, beta0, beta1)
end

maxnorm(x) = maximum(abs, x)



@testset "SCTransform.jl" begin
	@testset "500_PBMC_50genes" begin
		filename = "data/500_PBMC_50genes/filtered_feature_bc_matrix.h5"
		X,f,b = read10x(filename)

		pnonreg_ans = read_model_pars("data/500_PBMC_50genes/model_pars.csv")
		gm_ans = readdlm("data/500_PBMC_50genes/genemean.csv", Float64)

		pnonreg_ans2 = read_model_pars("data/500_PBMC_50genes/model_pars2.csv")
		bw_ans2 = 0.243867597143564
		preg_ans2 = read_model_pars("data/500_PBMC_50genes/model_pars_fit2.csv")
		gm_ans2 = readdlm("data/500_PBMC_50genes/genemean2.csv", Float64)
		Z_ans2 = readdlm("data/500_PBMC_50genes/transformed2_every_10th_cell.csv", ',', Float64)



		@testset "TableType=$TableType" for TableType in (NamedTuple, DataFrame)
			@testset "Non-regularized" begin
				pnonreg = scparams_estimate(TableType, X; chunk_size=20, feature_names=f.name)
				@test pnonreg isa TableType

				@test f.name[pnonreg.featureInd] == pnonreg_ans.genename

				@test pnonreg.logGeneMean ≈ gm_ans
				@test pnonreg.beta0_estimate ≈ pnonreg_ans.beta0
				@test pnonreg.beta1_estimate ≈ pnonreg_ans.beta1

				# The ground truth is really poor for theta, in particular for cases when the ml
				# estimate of theta is very large.
				# We thus test different genes at different levels of strictness.
				theta_ind_strict = setdiff(1:38, [3,7,10,12,13,15,16,19,20,31,32,36,38])
				theta_ind_approx = setdiff(1:38, [3,10,13,16,19,31,32,36])
				@test pnonreg.theta_estimate[theta_ind_strict] ≈ pnonreg_ans.theta[theta_ind_strict] rtol=0.01
				@test pnonreg.theta_estimate[theta_ind_approx] ≈ pnonreg_ans.theta[theta_ind_approx] atol=.5
				@test all(>(-0.01), pnonreg.theta_estimate .- pnonreg_ans.theta) # Our estimates should always be larger
			end

			# Testing regularization is difficult when we flag thetas as outliers differently.
			# We thus skip thetas flagged as outliers for this part.
			@testset "Regularized" begin
				f_ind = setdiff(1:50, [7,15,19,22,25,42,43,47])
				X2 = X[f_ind,:]
				f2 = TableType(SCTransform.subset_rows(f, f_ind)) # TODO: avoid using internal function

				pnonreg = scparams_estimate(TableType, X2; chunk_size=19, feature_names=f2.name)
				@test pnonreg isa TableType

				@test f.name[f_ind][pnonreg.featureInd] == pnonreg_ans2.genename

				@test pnonreg.logGeneMean ≈ gm_ans2
				@test pnonreg.beta0_estimate ≈ pnonreg_ans2.beta0
				@test pnonreg.beta1_estimate ≈ pnonreg_ans2.beta1
				@test pnonreg.theta_estimate ≈ pnonreg_ans2.theta atol=0.5


				pnonreg2 = scparams_detect_outliers(pnonreg)
				@test pnonreg2 isa TableType

				bw = scparams_bandwidth(pnonreg2)
				@test bw ≈ bw_ans2 rtol=0.1

				# --- Inexact comparisons below, since the bandwidths etc. differ a bit ---
				preg = scparams_regularize(pnonreg2,bw)
				@test preg isa TableType

				@test f.name[f_ind][preg.featureInd] == preg_ans2.genename

				@test preg.logGeneMean ≈ gm_ans2
				@test preg.beta0 ≈ preg_ans2.beta0 norm=maxnorm atol=0.2 # We might be able to tighten this
				@test preg.beta1 ≈ preg_ans2.beta1 norm=maxnorm atol=0.2 # We might be able to tighten this
				@test preg.theta ≈ preg_ans2.theta atol=0.5

				preg2 = scparams(X2, TableType(f2); use_cache=false) # whole pipeline in one call
				@test preg2 isa TableType

				@test preg2.name == pnonreg_ans2.genename
				@test !hasproperty(preg2,:featureInd)

				@test preg2.logGeneMean ≈ gm_ans2
				@test preg2.beta0 ≈ preg_ans2.beta0 norm=maxnorm atol=0.2 # We might be able to tighten this
				@test preg2.beta1 ≈ preg_ans2.beta1 norm=maxnorm atol=0.2 # We might be able to tighten this
				@test preg2.theta ≈ preg_ans2.theta atol=0.5

				Z2 = sctransform(X2, f2, preg2)
				@test Z2[:,1:10:end] ≈ Z_ans2 rtol=0.1
			end
		end

		@testset "sctransform" begin
			f_ind = setdiff(1:50, [7,15,19,22,25,42,43,47])
			X2 = X[f_ind,:]
			f2 = SCTransform.subset_rows(f, f_ind) # TODO: avoid using internal function

			row_ind = indexin(preg_ans2.genename,f.name)
			preg2 = (; id = f.id[row_ind],
			           feature_type = f.feature_type[row_ind],
			           beta0 = preg_ans2.beta0,
			           beta1 = preg_ans2.beta1,
			           theta = preg_ans2.theta)

			Z2 = sctransform(X2, f2, preg2)

			@test Z2[:,1:10:end] ≈ Z_ans2 rtol=1e-4
			@test Z2[:,1:10:end] ≈ Z_ans2 norm=maxnorm rtol=1e-4


			f_subset = 2:3:length(preg2.id)
			c_subset = 7:18:size(X2,2)

			Z2s = sctransform(X2, f2, SCTransform.subset_rows(preg2,f_subset))
			@test Z2s == Z2[f_subset,:]

			Z2s = sctransform(X2, f2, preg2; cell_ind=c_subset)
			@test Z2s == Z2[:,c_subset]

			Z2s = sctransform(X2, f2, SCTransform.subset_rows(preg2,f_subset); cell_ind=c_subset)
			@test Z2s == Z2[f_subset,c_subset]

			preg3 = (; id=[preg2.id[1], "NOT_A_GENE_ID"],
			           feature_type = preg2.feature_type[1:2],
			           beta0 = preg2.beta0[1:2],
			           beta1 = preg2.beta1[1:2],
			           theta = preg2.theta[1:2])
			@test_throws DomainError sctransform(X2, f2, preg3)
		end

	end
end

