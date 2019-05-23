
__all__ = ['isTs',]

def isTs(pair):
    if pair in ['AG','GA','CT','TC']:
    	return 1
    elif pair in ['AC','CA','AT','TA','GC','CG','GT','TG']:
    	return 0
    else: return None
