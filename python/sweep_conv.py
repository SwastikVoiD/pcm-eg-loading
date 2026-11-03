from sweep import run
for p in ['erythritol','paraffin']:
    for k in ['low','mid','high']:
        run(p,k,1000.0,conv=True)
print('CONV DONE',flush=True)
