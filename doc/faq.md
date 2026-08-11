<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!---->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# FAQ

## Snakemake `retrieve_*` rules fail with `SSL: CERTIFICATE_VERIFY_FAILED` error

If your workflow fails while downloading input data (for example during `snakemake -n`), and you may see errors like this:

```console
$ pixi run snakemake -n --cores 1
Building DAG of jobs...
ConnectError while get'ing https://data.pypsa.org/workflows/eur/corine/v18_5/manifest.yaml
[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: self-signed certificate in certificate chain (_ssl.c:1032)
...
httpx.ConnectError: [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: self-signed certificate in certificate chain (_ssl.c:1032)
```

This typically happens on networks where a university/company proxy intercepts HTTPS traffic and presents a different and incompletely configured SSL certificate.

To inspect the certificate chain you actually receive for `data.pypsa.org`, run:

```console
pixi run openssl s_client -connect data.pypsa.org:443 -servername data.pypsa.org -showcerts </dev/null 2>/dev/null \
| awk '/BEGIN CERTIFICATE/{i++} i{print > "/tmp/chain-cert-" i ".pem"} END{for (j=1;j<=i;j++) print "/tmp/chain-cert-" j ".pem"}' \
| while read -r cert; do
        echo "=== $cert ==="
        openssl x509 -in "$cert" -noout -subject -issuer
        openssl x509 -in "$cert" -noout -text \
            | awk '
                /Version:/ {print}
                /Authority Key Identifier/ {print; getline; print}
                /Subject Key Identifier/ {print; getline; print}
            '
    done
```

This certificate should be for `data.pypsa.org` and issued by a trusted certificate authority.
If it is not and e.g. self-signed or issued by your organisation, then this can cause the error above.

If the certificate chain is replaced by a proxy certificate, common fixes are:

1. Run PyPSA-Eur from a network without that proxy (for example from home).
2. Ask your IT department for help, they can either provide a proxy bypass or ensure that their proxy certificate is correctly configured and trusted by your system.

This is usually faster than trying to disable certificate verification or hot fixes.