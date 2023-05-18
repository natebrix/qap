import pandas as pdAAA
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn import tree
from sklearn.model_selection import train_test_split
import lightgbm as lgb
from sklearn.metrics import mean_absolute_percentage_error
import datetime

# This might be a nice "intermediate strategy". That is:
#
#  Rule 4 is best at the top
#  Rule 2 is best at the bottom (when we are starting to fathom a lot)
#  ---> Rule 5 is best in the middle, even though it is also (as implemented) looking at fathoms.

# Let's run some size 20 problems down to a certain depth, starting with a size 25 root.


# So, what does "the middle" mean?
#   It's when we think there is going to be a 'significant' percentage of children fathomed on this level.
#   Let X be that percentage. Then that means:
#       b + u  >= z
#   where b is the bound, u is an estimate of the u value *AT THIS LEVEL* and z is the incument
#
#  or I can just simplify and use the gap ratio as I normally do.

# add a check for "B" matrix

# i need a reliable way to log B

# back to basics
# non-ML model that is better than rule 2.
#  - read cofficients from file?
#  - do a grid search and mess with functional form
#
# - I do see some very modest improvements by reading coefficients from a file and weighting
#    'the ends' more.
#   I see it when I replace rule 2 with rule 5 in nug20.par. (so lower down in the tree)
# how much mileage can I get out of just using the values of U anyway?!?
#
# basic problem once more:
#  - it seems to be hard to do much better than rule 2 when we are far down in the tree.
#        - how do I know this? What is the actual opportunity to do better?
#  - it seems to be hard to do better than rule 4
#        - how do I know this?
#  - this may help only for really big problems?!?

# another completely different idea:
#  - switch branching rule based on the percentage of children fathomed at previous level

# I can just compute arrays of things and do the spearman coefficient
# ss.spearmanr()

# If I really want to have a good result, I need to do better than rule 4.
# which means I need to log the results of rule 4.
#
# A good start would be to log B.

# maybe u variance is better!!!
#
#>>> xv=data2['x_i'].apply(np.var)
#>>> ss.spearmanr(xv.argsort().argsort(), np.array(data2['node'].argsort().argsort()))
#SpearmanrResult(correlation=0.059340659340659345, pvalue=0.8403017626906826)
#>>> ss.spearmanr(uv.argsort().argsort(), np.array(data2['node'].argsort().argsort()))
#SpearmanrResult(correlation=0.2351648351648352, pvalue=0.4183338961985912)
#>>> ss.spearmanr(us.argsort().argsort(), np.array(data2['node'].argsort().argsort()))
#SpearmanrResult(correlation=0.14285714285714285, pvalue=0.6261174762253241)

def rank(v):
    return np.array(v.argsort().argsort())

def spearman(data):
    xv = data['x_i'].apply(np.var)
    us = data['u_i'].apply(sum)
    uv = data['u_i'].apply(np.var)
    best = data['node']
    s_us = ss.spearmanr(us, best)
    s_uv = ss.spearmanr(uv, best)
    s_xv = ss.spearmanr(xv, best)
    print(f'x variance: {s_xv}')
    print(f'u variance: {s_uv}')
    print(f'u sum     : {s_us}')

# another thing we could do is try to learn B.
# that is, whatever score we would get from lookahead branching?!?

# verification steps:
#  1. feature validation
#  2. label validation
#       - improved label and eye validated
#  3. model parameter validation
#       - have not done this
#  4. score validation
#  5. offline stats validation
#  6. code validation
#  7. integration validation

# if a child is eliminated, it's value is very large...
# could think about percentage gap reduction (maxing at 100%)
#  
# min (1, u / (inc - bound - u))

# theories:
#  - it will be hard to outperform strong branching unless we add more features
#  - if we cannot outperform strong branching then this is a good rule for the middle/bottom of tree
#  - it might be hard to outperform rule 2 for the bottom of the tree
#  - therefore it might be best for the middle of the tree

# need to see if I can get better results.
#  -- need to test more in python
#  -- nned to verify the inputs and outputs precisely

# I can diversify my training data now.
#
# If I can learn better than lookahead branching, I will be in business.
#   -> maybe I can try to learn what is going on k levels down??
#   -> how many nodes do I have on level K?

# Any other feature ideas????
#    x?
#    bound gap?
#    C dominance?

# I guess I need to set the label gain.
# This tells me how much I should care about 
# I want the 'rank' to be high for 'good' things

# NO I think it is good to have lower rank = better. so rank 1 is always best

#>>> lg = [pow(2, i) for i in range(34)]
#>>> estimate(data, params={'label_gain':lg})

lg =  [pow(2, i) for i in range(34)] # consider changing


# stretch the values in a_small into an array of size n_big.
def stretch(a_small, n_big):
    n_small = len(a_small)
    ii = [int(i * n_small / n_big) for i in range(n_big)]
    return a_small[ii]

def squish(a_big, n_small):
    return np.concatenate([a_big[:n_small-n_small//2], a_big[-n_small//2:]])

def force_to_size(a, n):
    return squish(a, n) if len(a) > n else stretch(a, n)


# train_data = lgb.Dataset(data, label=label, feature_name=['c1', 'c2', 'c3'], categorical_feature=['c3'])

# train_data = lgb.Dataset(data, label=label, weight=w)
# param['metric'] = ['auc', 'binary_logloss']


def read_log(filename):
    df = pd.read_csv(filename, delimiter="|", header=None)
    df2 = df.rename(columns={2:'trial', 1:'rc', 3:'index', 0:'op'})
    df2['trial'] = df2['trial'].astype(int)
    return df2

def make_data_from_logs(logs):
    dfs = [read_log(log) for log in logs]
    datas = [make_data(df) for df in dfs]
    trial_offsets = np.concatenate(([0], np.cumsum([data['trial'].max() for data in datas])))
    for i, data in enumerate(datas):
        data['trial'] += trial_offsets[i]
    return pd.concat(datas)    

# more things to try:
#  lightgbm
#  reduce the tree size
#  observation weights
#  custom loss function

# r is a row.
def row_to_mat(r):
    n = int(np.sqrt(r.array.shape[0]))
    return np.reshape(np.array(r.loc[:]).astype(float), [n, n])


def reshape_matrix(u, name):
    ur = pd.DataFrame(index=u.index)
    ur[f'{name.lower()}'] = u.apply(row_to_mat, axis=1)
    return ur

def get_matrix_row_col(r, name='u'):
    u = r[name]
    rc = r['rc']
    index = int(r['index'])
    n = u.shape[0]
    return u[index, :] if rc == 'r' else u[:, index]


def make_matrix(df2, name):
    n = int(np.sqrt(df2.shape[1]-3))

    # get U entries and columns
    du = df2[(df2['op'] == name)]
    du = du[['trial'] + [i for i in range(4,4+n*n)]]
    du = du.rename(columns={i : name+str(i-4) for i in range(4,4+n*n)})
    du1 = du.set_index('trial')
    return reshape_matrix(du1, name)

def make_r(df2):
    dr = df2[(df2['op']=='t')]
    dr = dr.rename(columns={5:'inc', 6:'bound', 7:'time',  8:'node', 9:'fathom', 10:'enumerate'})
    dr['index'] = dr['index'].astype(int)
    dr = dr[['op','rc','trial','index','inc','bound','time','node','fathom','enumerate',11]]
    dr1 = dr.set_index('trial')
    return dr1

def gap_reduction(data):
    return data.apply(lambda r: np.minimum(1, r['u_raw'] / (r['inc'] - r['bound'])), axis=1)

def mean_center(data):
    raw = 'u_raw'
    return data.apply(lambda r: r[raw] / r[raw].mean(), axis=1)
    
def var_names(prefix, n):
    return [f'{prefix}{i}' for i in range(n)]

def array_col_to_df(data, col, prefix, n_features):
    return pd.DataFrame(data=unpack_array(data, col), columns=var_names(prefix, n_features), index=data.index)

def fathom_count(data):
    return data.apply(lambda r: (r['u_raw'] + r['bound'] >= r['inc']).sum(), axis=1)

def col_rank(data, score, ascending=False):
    return data.groupby(['trial'])[score].rank(ascending=ascending, method='first').astype(int)

def col_matrix_sum(data, name, n_features):
    data[f'{name}_i_raw'] = data.apply(lambda r: get_matrix_row_col(r, name=name), axis=1)
    data[f'{name}_i'] = data[f'{name}_i_raw'].apply(lambda r: force_to_size(r, n_features)).apply(np.sort)
    data[f'{name}_sum'] = data[f'{name}_i'].apply(sum)
    return data

def col_sum_scores(data, col_sum):
    data = data.join(data.groupby(['trial'])[f'{col_sum}'].max(), on='trial', rsuffix="_best")
    # u_sum_score is ratio between row sum and best for trial.
    data[f'{col_sum}_score'] = data[f'{col_sum}'] / data[f'{col_sum}_best']
    data[f'{col_sum}_score_rank'] = col_rank(data, f'{col_sum}_score', ascending=True)
    return data

def make_data(df, n_features=16, y_col='node'):
    r = make_r(df)
    u = make_matrix(df, 'U')
    x = make_matrix(df, 'X')
    b = make_matrix(df, 'B')

    data = r.join(u).join(x)
    if b.shape[0] > 0:
        data = data.join(b)
        data = col_matrix_sum(data, 'b', n_features)
        data = col_sum_scores(data, 'b_sum')
    data['u_raw'] = data['u'].copy()
    data['n_fathom'] = fathom_count(data)
    # This is the ratio by which u_ij reduces the gap.
    #    25 / (1150 - 1100 - 25) = 25 / 25) = 1 ==> wrong
    #data['u'] = gap_reduction(data) # todo no good: just a scaling
    data['u'] = mean_center(data)
    data = data.reset_index()
    data['choice'] = np.where(data['rc']=='r', data['index'].astype(int), 10+data['index'].astype(int))
    
    data['u_i_raw'] = data.apply(lambda r: get_matrix_row_col(r, name='u'), axis=1)
    data['u_i'] = data['u_i_raw'].apply(lambda r: force_to_size(r, n_features)).apply(np.sort)
    data['u_sum'] = data['u_i'].apply(sum)
    data_u = array_col_to_df(data, 'u_i', 'u', n_features)
    #data_u = pd.DataFrame(data=unpack_array(data, 'u_i'), columns=var_names('u', n_features), index=data.index)

    data['x_i_raw'] = data.apply(lambda r: get_matrix_row_col(r, name='x'), axis=1)
    data['x_i'] = data['x_i_raw'].apply(lambda r: force_to_size(r, n_features)).apply(np.sort)
    data_x = array_col_to_df(data, 'x_i', 'x', n_features)

    d_list = [data, data_u, data_x]
    if b.shape[0] > 0:
        d_list += [array_col_to_df(data, 'b_i', 'b', n_features)]
    data = pd.concat(d_list, axis=1)

    # get the best row/col sums
    data = data.join(data.groupby(['trial'])['u_sum'].max(), on='trial', rsuffix="_best")
    # u_sum_score is ratio between row sum and best for trial.
    data['u_sum_score'] = data['u_sum'] / data['u_sum_best']
    data['u_sum_score_rank'] = col_rank(data, 'u_sum_score', ascending=True)

    # suppose min is 500 and I am 750. Then y = 500/750 = .6666
    # alternatively could do (x - min) / min. Then y = (750 - 500) / 500 = 0.5
    # altneratively could do (x - min) / (max - min) and have labels from 0 - 1
    # data['y'] = data.groupby('trial')[y_col].transform('min') / data[y_col]
    by_trial_y = data.groupby('trial')[y_col]
    data['y_min'] = by_trial_y.transform('min')
    data['y_max'] = by_trial_y.transform('max')
    data['y'] = (data[y_col] - data['y_min']) / (data['y_max'] - data['y_min']) 
    data['y0'] = data[y_col] / data['y_min']
    data['y_rank'] = col_rank(data, 'y', ascending=True) # higher rank means better
    return data

def unpack_array(data, col):
    return np.stack(data.loc[:, col].array)


def node_count(data, rank_col, op=max):
    idx = data.groupby(['trial'])[rank_col].transform(op) == data[rank_col]
    nc = data[idx]['node']
    #nc = data[data[rank_col]==1]['node']  # wrong beause I want the max in each group
    trials = data['trial'].unique()
    if nc.shape[0] != len(trials):
        raise ValueError('Incorrect number of records in node count')
    return nc.sum()


#######

def lgb_train(data_train, X_train, y_train, params):
    node_train = data_train['node'] # slightly better with weights
    train_data = lgb.Dataset(X_train, label=y_train, weight=node_train)
    #train_data = lgb.Dataset(X_train, label=y_train)
    bst = lgb.train(params, train_data)
    return bst

def tree_train(data_train, X_train, y_train, params):
    return tree.DecisionTreeRegressor().fit(X_train, y_train)

# larger rank means better.
def lgb_rank_train(data_train, X_train, y_train, params):
    # todo class_weight
    gids = data_train.groupby("trial")["trial"].count().to_numpy() # size of each group
    rnk = lgb.LGBMRanker(**params)
    rnk.fit(X_train, y_train, group=gids)
    return rnk

#######


def model_stats(data, base_score, scores):
    n_b = node_count(data, base_score)
    n = [node_count(data, name) for name in scores]
    by_trial = data.groupby(['trial'])['node']
    n_best = by_trial.min().sum()
    n_worst = by_trial.max().sum()
    n_med = by_trial.median().sum()
    print('---------------------------------')
    print('           MODEL STATS           ')
    print('---------------------------------')
    print(f'size      \t\t= {data.shape}')
    print(f'worst     \t\t= {n_worst:8.1f}')
    print(f'median    \t\t= {n_med:8.1f}')
    print(f'{base_score}\t= {n_b:8.1f}')
    for i, name in enumerate(scores):
        n_u = node_count(data, name)
        ratio = (n_u - n_b) / n_u # improvement from baseline
        print(f'{name}\t= {n_u:8.1f}')
        print(f' ratio     \t\t= {ratio:8.3f}')
    print(f'best           \t\t= {n_best:8.1f}')

    
def feature_cols(n_features):
    return var_names('u', n_features) # + ['u_sum_score'] + var_names('x', n_features)

# Split a data set into train/test, keeping all rows from the same trial together
def split_by_trial(data, test_size):
    trials = data['trial'].unique()
    np.random.shuffle(trials)
    cutoff = len(trials) - int(len(trials) * test_size)
    
    data_train = data[data['trial'].isin(trials[:cutoff])].copy()
    data_test = data[data['trial'].isin(trials[cutoff:])].copy()
    return data_train, data_test

# get all of the rank columns in a data set, except those in the exclude list.
def rank_columns(data, exclude):
    return [c for c in data.columns if ('score_rank' in str(c)) and (str(c) not in exclude)]

# evaluate a model to produce row/column rankings for branching. The provided data will be split
# into test and training sets. If no model is provided, then one will be trained using the specified
# train function, passing any supplied parameters (in params).
#
# Either way, once a trained model is obtained, it will be used to predict ranks for the test set.
# Then the quality of the ranks will be evaluated and compared to a 'baseline' ranking.
def estimate(data, y_col='y_rank', train=lgb_rank_train, params={}, model=None, test_size = 0.2):
    print(f'estiamte: {train.__name__} -> {y_col} with {params}')
    n = data.loc[0, 'u_i'].shape[0]
    if test_size < 1.0:
        print(f'Splitting into train, test with ratio {test_size}')
        data_train, data_test = split_by_trial(data, test_size=test_size)
    else:
        data_train = data
        data_test = data

    X_train = data_train[feature_cols(n)]
    y_train = data_train[y_col]
    X_test = data_test[feature_cols(n)]
    y_test = data_test[y_col]

    if model is None:
        print('Training model')
        model = train(data_train, X_train, y_train, params)
    y_pred_all = model.predict(data[feature_cols(n)])
    y_pred = model.predict(X_test)
    
    e_sum = mean_absolute_percentage_error(y_test, data_test['u_sum_score'])
    e_mod = mean_absolute_percentage_error(y_test, y_pred)
    print(f'u_sum_mape = {e_sum:.3f}, model_mape = {e_mod:.3f}')
    data['model_score'] = y_pred_all
    data['model_score_rank'] = col_rank(data, 'model_score', ascending=True)
    # todo I want this to work for both regular prediction and ranking.
    # for regular prediction, higher is better.
    # for ranking, lower is better.
    # I could fix this by having a wrapper class for the ranker, so when we call predict
    # it returns the negated.
    data_test['model_score'] = y_pred # todo this is inconsistent between std and rank
    data_test['model_score_rank'] = col_rank(data_test, 'model_score', ascending=True)
    b_s = 'u_sum_score_rank'
    model_stats(data_test, base_score=b_s, scores=rank_columns(data_test, exclude=[b_s]))
    return model, data_test


def estimate_from_file(filename, y_col='y_rank', train=lgb_rank_train, params={}, model=None):
    df = read_log(filename)
    data = make_data(df)
    return estimate(data, y_col=y_col, train=train, params=params, model=model)


class LGBGenerateC:
    def __init__(self, path, indent=0):
        self.path = path
        self.indent = indent
     
    def __enter__(self):
        self.file = open(self.path, 'w')
        return self
 
    def __exit__(self, *args):
        self.file.close()

    def write(text):
        pass

    def arg_list(self, args):
        t = [f'{a["type"]} {a["name"]}' for a in args]
        return ','.join(t)

    def comment(self, text):
        self.file.write( f'{" "*self.indent}// {text}\n')
        
    def function_start(self, name, args, return_type='void'):
        a_list = self.arg_list(args)
        self.file.write( f'{" "*self.indent}{return_type} {name}({a_list}) {{\n')
        self.indent += 2

    def block_end(self):
        self.indent -= 2
        self.file.write(f'{" "*self.indent}}}\n')
    
    def function_end(self):
        self.block_end()
        self.file.write('\n')
        
    def if_start(self, clause):
        self.file.write(f'{" "*self.indent}if ({clause}) {{\n')
        self.indent += 2

    def if_end(self):
        self.indent -=2
        self.file.write(f'{" "*self.indent}}}\n')

    def else_start(self):
        self.file.write(f'{" "*self.indent}else {{\n')
        self.indent += 2

    def else_end(self):
        self.indent -= 2
        self.file.write(f'{" "*self.indent}}}\n')

    def return_clause(self, clause):
        self.file.write(f'{" "*self.indent}return {clause};\n')

    def sum_expr(self, items):
        sum_items = ' + '.join(items)
        return f'{sum_items}'
        
class LGBGeneratePython:
    def __init__(self, path, indent=0):
        self.path = path
        self.indent = indent
     
    def __enter__(self):
        self.file = open(self.path, 'w')
        return self
 
    def __exit__(self, *args):
        self.file.close()

    def write(text):
        pass

    def arg_list(self, args):
        t = [f'{a["name"]}' for a in args]
        return ','.join(t)

    def comment(self, text):
        self.file.write( f'{" "*self.indent}# {text}\n')

    def function_start(self, name, args, return_type='void'):
        a_list = self.arg_list(args)
        self.file.write( f'{" "*self.indent}def {name}({a_list}):\n')
        self.indent += 2

    def block_end(self):
        self.indent -= 2
        self.file.write(f'\n')
    
    def function_end(self):
        self.block_end()
        self.file.write('\n')
        
    def if_start(self, clause):
        self.file.write(f'{" "*self.indent}if {clause}: \n')
        self.indent += 2

    def if_end(self):
        self.indent -=2

    def else_start(self):
        self.file.write(f'{" "*self.indent}else:\n')
        self.indent += 2

    def else_end(self):
        self.indent -= 2

    def return_clause(self, clause):
        self.file.write(f'{" "*self.indent}return {clause}\n')

    def sum_expr(self, items):
        sum_expr = '[' + ' , '.join(items) + ']'
        return f'sum({sum_expr})'
        
def print_tree(t, gen):
    if 'split_feature' in t:
        gen.if_start(f'u[{t["split_feature"]}] <= {t["threshold"]}')
        if 'left_child' in t:
            print_tree(t['left_child'], gen)
        gen.if_end()
        if 'right_child' in t:
            gen.else_start()
            print_tree(t['right_child'], gen)
            gen.else_end()
    else:
        gen.return_clause(f'{t["leaf_value"]}')

def get_generate(language, filename):
    return LGBGenerateC(filename) if language.lower()=='c' else LGBGeneratePython(filename)
        
def model_to_code(model, name='lgb_train', prefix='u', language='c'):
    t = model.dump_model() if hasattr(model, 'dump_model') else model.booster_.dump_model()
    filename = f'{name}.{language}'
    print(f'writing model to {filename}, language = {language}')
    with get_generate(language, filename) as gen:
        f_args = [{'type':'double *', 'name':f'{prefix}'}]
        fs = []
        gen.comment(f'Generated on {datetime.datetime.now()}');
        for i, b in enumerate(t['tree_info']):
            fs += [f'boost_{i}({prefix})']
            gen.function_start(name=f'boost_{i}', args=f_args, return_type='double')
            print_tree(b['tree_structure'], gen)
            gen.function_end()
        gen.function_start(name='boost', return_type='double', args=f_args)
        gen.return_clause(gen.sum_expr(fs))
        gen.function_end()
            

def predict_str(model, s):
    d = np.array(list(map(float, ss.split(' '))))
    return model.predict([d])


# this is what we would have chosen at the current node.
# what we want to know is where this ranks within the y's.
#u2['sum_c'] = u2.apply(lambda r: r['U'].sum(axis=0), axis=1)
#u2['sum_r'] = u2.apply(lambda r: r['U'].sum(axis=1), axis=1)
#u2['i_c'] = u2.apply(lambda r: np.argmax(r['sum_c']), axis=1)
#u2['i_r'] = u2.apply(lambda r: np.argmax(r['sum_r']), axis=1)
