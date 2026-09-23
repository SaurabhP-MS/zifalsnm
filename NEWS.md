# zifalsnm 0.1.0

* Updated implementation of the CAVI algorithm for the ZIFA-LSNM model, used for the revised manuscript.
* The last column of the count matrix is now used as the reference taxon; `Estimated_Compositions` contains all `P` taxa.
* New arguments to `ZIFA_LSNM()`: `number_of_factors` (replaces `num_fac`), `epsilon`, `NU1Prior`, `NU2Prior`, `G1Prior`, `G2Prior`, `AlphaPrior` (replace `NU1P`, `NU2P`, `G1P`, `G2P`, `AlphaP`), `Max_Iter` and `verbose`.
* Convergence is now based on the relative change in the ELBO.
* The output now includes the trace of every variational parameter across iterations.
* Updated reproducibility scripts in `Reproducible_R_Codes/`.

# zifalsnm 0.0.0.9000

* First version.
