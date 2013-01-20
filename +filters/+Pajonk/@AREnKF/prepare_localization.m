function [] = prepare_localization(this, model, H)

dist = model.distanceMatrix();
dl = model.decorrelationLength();

this.covLoc = this.opts.localization_function(dist, dl*this.opts.localization_function_width, this.N);
this.mCovLoc = H*this.covLoc;

end

