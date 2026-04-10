# distutils: language = c++
# cython: language_level=3

from cpython cimport array
import array

from nlbayes.ModelORNOR cimport *
import pandas as pd
import numpy as np
from tqdm import tqdm
from cpython.exc cimport PyErr_CheckSignals
from pprint import pprint


cdef class PyModelORNOR:
    cdef ModelORNOR *c_model
    cdef public object _seeds

    def __cinit__(self,
                  network,
                  evidence=dict(),
                  set active_tf_set=set(),
                  uniform_t=False,
                  t_alpha=None,
                  t_beta=None,
                  zy=0.,
                  zn=0.,
                  s_leniency=0.1,
                  n_graphs=3,
                  verbosity=0,
                  base_seed=5489):

        # ---- params / defaults ----
        zy_value = 0.99 if zy == 0. else float(zy)

        # Compute zn if evidence provided and user didn't set it
        if len(evidence) > 0 and zn == 0.:
            # IMPORTANT: network is a dict-of-dicts, so values() is required
            n_edges = sum(len(d) for d in network.values())
            n_edges_deg = 0
            for src, trg_dict in network.items():
                for trg, mor in trg_dict.items():
                    if (trg in evidence.keys() and evidence[trg] != 0):
                        n_edges_deg += 1
            zn_value = (n_edges_deg / n_edges / 10) if n_edges > 0 else 1.0
        elif zn == 0.:
            zn_value = 1.
        else:
            zn_value = float(zn)

        # Theta prior
        if len(evidence) == 0:
            if t_alpha is None and t_beta is None:
                t_alpha = 18.
                t_beta = 2.
            else:
                assert t_alpha is not None and t_beta is not None, (
                    "Must provide both t_alpha and t_beta.")
        elif uniform_t:
            t_alpha = 1.
            t_beta = 1.
        else:
            if t_alpha is None and t_beta is None:
                t_alpha = 2.
                t_beta = 2.
            else:
                assert t_alpha is not None and t_beta is not None, (
                    "Must provide both t_alpha and t_beta.")

        if verbosity >= 2:
            print("\tT alpha   : ", t_alpha)
            print("\tT beta    : ", t_beta)
            print("\tS leniency: ", s_leniency)
            print("\tZY value  : ", zy_value)
            print("\tZN value  : ", zn_value)
            print("\t# Graphs  : ", n_graphs)
            print("\tbase_seed : ", base_seed)

        # ---- build network_c (vector of edges) ----
        cdef src_trg_pair_t src_trg
        cdef network_edge_t network_edge
        cdef network_t network_c = network_t()

        for src_uid, trg_dict in network.items():
            for trg_uid, mor in trg_dict.items():
                src_trg = (str(src_uid).encode("utf8"),
                           str(trg_uid).encode("utf8"))
                network_edge = src_trg, int(mor)
                network_c.push_back(network_edge)

        # Determine which genes are present in the regulatory network
        # IMPORTANT: values() is required (values is a method on dict)
        genes_in_network = set([g for r in network.values() for g in r.keys()])

        # ---- build evidence_c (map) ----
        cdef evidence_dict_t evidence_c = evidence_dict_t()
        for trg_uid, trg_de in evidence.items():
            if trg_uid in genes_in_network:
                evidence_c.insert((str(trg_uid).encode("utf8"), int(trg_de)))

        # ---- build active_tf_set_c (set) ----
        cdef prior_active_tf_set_t active_tf_set_c = prior_active_tf_set_t()
        for src_uid in active_tf_set:
            active_tf_set_c.insert(str(src_uid).encode("utf8"))

        # ---- sprior ----
        sprior = [
            1.0 - s_leniency, 0.9 * s_leniency, 0.1 * s_leniency,
            0.5 * s_leniency, 1.0 - s_leniency, 0.5 * s_leniency,
            0.1 * s_leniency, 0.9 * s_leniency, 1.0 - s_leniency
        ]
        cdef double[9] sprior_c = array.array("d", sprior)

        # ---- construct C++ model ----
        self.c_model = new ModelORNOR(
            network_c, evidence_c, active_tf_set_c,
            sprior_c,
            float(t_alpha), float(t_beta),
            float(zy_value), float(zn_value),
            <unsigned int> n_graphs,
            <unsigned int> base_seed
        )

        # Cache seeds so adapter can read them safely
        try:
            self._seeds = [int(x) for x in self.c_model.get_seeds()]
        except Exception:
            self._seeds = []


    def __dealloc__(self):
        if self.c_model is not NULL:
            del self.c_model
            self.c_model = NULL


    # -----------------
    # Convergence helpers
    # -----------------
    def get_gelman_rubin(self):
        gr_list = self.c_model.get_gelman_rubin()
        return [id.decode("utf8").split("_")[:2] + [gr] for id, gr in gr_list]

    def get_max_gelman_rubin(self):
        return self.c_model.get_max_gelman_rubin()

    @property
    def seeds(self):
        return list(self._seeds)


    # -----------------
    # Sampling
    # -----------------
    def sample_posterior(self, N, gr_level, burnin=False, show_progress=True):
        if burnin:
            print("\nInitializing model burn-in ...")
            converged = False
            while not converged:
                status = self.sample_n(200, 20, 5.0, show_progress)
                converged = status == 0
                if status == -1:
                    print("Interrupt signal received")
                    return
            self.burn_stats()
            print("Burn-in complete ...")

        n_sampled = 0
        converged = False
        i = 1
        while n_sampled < N and not converged:
            n = min(200 * i, N - n_sampled)
            status = self.sample_n(n, 5, gr_level, show_progress)
            converged = status == 0
            n_sampled += n
            i += 1
            if status == -1:
                print("Interrupt signal received")
                return

        if not converged:
            x = dict(self.c_model.get_gelman_rubin())
            n_vars = len(x)
            x = {k: v for k, v in x.items() if v > gr_level}
            n_did_not_converge = len(x)
            print("\nThere are", n_vars, "random variables in the model.",
                  n_did_not_converge, "of them did not converge.")
            if n_did_not_converge < 20:
                pprint(x)

    def sample_n(self, N, dN, gr_level, show_progress=True, quiet=False):
        if not quiet:
            print()
        gr = float("inf")
        status = 0
        n = 0

        try:
            if show_progress:
                progress = tqdm(total=N)
            while n < N and gr > gr_level:
                dN = min(dN, N - n)
                if show_progress:
                    progress.update(dN)

                self.c_model.sample_n(dN)
                n += dN
                gr = self.c_model.get_max_gelman_rubin()
                PyErr_CheckSignals()

            if show_progress:
                progress.total = n
                progress.update(0)

        except KeyboardInterrupt:
            status = -1
        finally:
            if show_progress:
                progress.close()

        converged = gr <= gr_level

        if converged:
            if not quiet:
                print("Converged after", self.c_model.total_sampled, "samples")
        elif status == 0:
            status = 1
            if not quiet:
                print("Drawed", self.c_model.total_sampled, "samples so far")
        elif status == -1:
            if not quiet:
                print("\nProcess interrupted.")
                print("Drawed", self.c_model.total_sampled, "samples so far")

        if not quiet:
            print("Max Gelman-Rubin statistic is", gr, "(target",
                  "was" if converged else "is", gr_level, ")")

        return status

    def burn_stats(self):
        self.c_model.burn_stats()


    # -----------------
    # Posterior extraction (what your adapter expects)
    # -----------------
    def inference_posterior_df(self, annotation=dict()):
        """
        Returns a DataFrame with TF posterior mean/sd.
        annotation: dict like { "TF_ID": "TF_NAME" } (optional)
        """
        means = self.c_model.get_posterior_means("S")
        sdevs = dict(self.c_model.get_posterior_sdevs("S"))

        rows = []
        for var_id, mu in means:
            vid = var_id.decode("utf8")
            tf_id = vid
            tf_name = annotation.get(tf_id, tf_id)
            sd = float(sdevs.get(var_id, np.nan))
            rows.append((tf_id, tf_name, float(mu), sd))

        df = pd.DataFrame(rows, columns=["tf_id", "tf", "mean", "sd"])
        df.sort_values(["mean"], ascending=False, inplace=True, ignore_index=True)
        return df


    # -----------------
    # Properties
    # -----------------
    @property
    def network(self):
        cdef int mor
        cdef std_string src_c, trg_c
        out = {}
        for i in range(self.c_model.network.size()):
            (src_c, trg_c), mor = self.c_model.network[i]
            src = src_c.decode()
            trg = trg_c.decode()
            if src not in out:
                out[src] = {}
            out[src][trg] = mor
        return out

    @property
    def evidence(self):
        return {uid.decode(): val for uid, val in self.c_model.evidence}

    @property
    def active_tf_set(self):
        return {uid.decode() for uid in self.c_model.active_tf_set}

    @property
    def _config(self):
        return {
            "t_alpha": self.c_model.t_alpha,
            "t_beta": self.c_model.t_beta,
            "zy": self.c_model.zy,
            "zn": self.c_model.zn,
            "n_graphs": self.c_model.n_graphs,
            "seeds": list(self._seeds),
        }
