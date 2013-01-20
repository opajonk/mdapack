function [ results ] = singleRun( parms )
%SINGLERUN Performs a single data assimilation experiment

if (parms.progress)
    fprintf(1, 'progress:');
end

method = parms.method();
model = parms.model();
representation = parms.representation(model);

stat = statistics.Statistics(parms, model, representation, method);

for t=1:model.timeStepCount()
    stat.storeForecast(t);
    
    if (model.hasMeasurement())
        method.update(model, representation);
        stat.storeFilterStatistics(t);
        
        if (parms.progress)
            fprintf(1, '.');
        end
    end
    
    stat.storeAnalysis(t);
    
    model.step(representation, t);
    
    if (parms.visualize)
        model.plot(representation);
        drawnow;
    end
end

if (parms.progress)
    fprintf(1,'done\n');
end

results = stat.getResults();
end