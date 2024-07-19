import os

class Error (Exception): pass
class Summary:
    def __init__(self,out_dir):
        self.out_dir = out_dir
    @staticmethod
    def summarise(out_dir):
        folders = sorted(os.listdir(out_dir))
        with open('summary.csv','w') as fobj:
            for d in folders:
                pred = os.path.join(out_dir,d,'pred.csv')
                if os.path.isfile(pred):
                    with open(pred, 'r') as f:
                        next(f)
                        fobj.write(next(f))

    def run(self):
        Summary.summarise(self.out_dir)
