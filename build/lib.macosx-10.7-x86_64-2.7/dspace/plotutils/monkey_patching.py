'''
'''

def monkeypatch_method(classes):
    '''
    '''
    def decorator(func):
        to_patch = classes
        if isinstance(to_patch, list) is False:
            to_patch = [classes]
        for cls in to_patch:
            setattr(cls, func.__name__, func)
        return func
    return decorator

